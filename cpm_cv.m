function Prediction = cpm_cv(Subjects_data,Subjects_score,pthresh,k_folds,corr_type,Model_Method)
% Runs cross validation for CPM
%   REQUIRED INPUTS
%        Subjects_data:        
%                MxMxN matrix containing all individual-subject connectivity matrices
%                m is number of nodes in the chosen brain atlas
%                n is the number of subjects
%------------------------------------------------------------------------------------%
%        Subjects_score:            
%               the continuous variable to be predicted for all subjects(e.g., behavioral scores) 
%               Allowed dimensions are 2D (i x num_subjects)
%------------------------------------------------------------------------------------%
%        pthresh:      
%               p-value threshold for feature selection,0.05 or 0.01 is common
%------------------------------------------------------------------------------------%
%        k_folds:    
%               Number of partitions for dividing the sample,5 or 10 is common
%               (e.g., 2 =split half, 10 = ten fold)
%               if kfolds == num_subjects % doing leave-one-out
%------------------------------------------------------------------------------------%
%        no_iterations:
%                 Repeatition times for the random n-folds cross-validation
%------------------------------------------------------------------------------------%
%        corr_type:
%                correlate all edges with behavior：     
%                    'Pearson' or 'Spearman' or 'Partial' or 'Robust'
%                'Pearson' 
%                    %correlate all edges with Subjects_score using Pearson correlation
%                'Spearman'
%                    %correlate all edges with Subjects_score using rank correlation
%                'Partial'
%                    %correlate all edges with Subjects_score using partial correlation
%                'Robust'
%                     %correlate all edges with Subjects_score using robust regression
%-------------------------------------------------------------------------------------% 
%          Mask_Method: 
%                 set threshold and define masks
%                 0 or 1; if 1, creating a weighted mask using sigmoidal function
%---------------------------------------------------------------------------------%
%          Model_Method: 
%                 build model on TRAIN subs and run model on TEST sub
%                 0 or 1; if 1, building model on TRAIN subs combining both postive
%                 and negative features  and run model on TEST sub
%---------------------------------------------------------------------------------%
num_subjects = size(Subjects_data,3); % get number of subjects 
number_node = size(Subjects_data,1);% get number of  nodes in the chosen brain atlas
RandID = randperm(num_subjects); 
EachFold_Quantity = floor(num_subjects/k_folds);
%---------------------------------------------=
fprintf('\n# Running over %1.0f Folds.\nPerforming fold no. ',k_folds);
test_score = [];
pred_test_score = [];
for leftout = 1:k_folds
    fprintf('%1.0f ',leftout);    
    if k_folds == num_subjects % doing leave-one-out
        testinds = RandID(leftout);
        traininds = setdiff(RandID,testinds);
    else
        Origin_ID{leftout} = RandID([(leftout - 1) * EachFold_Quantity + 1: leftout * EachFold_Quantity])';
        Reamin = mod(num_subjects, k_folds);
        for j = 1:Reamin
            Origin_ID{j} = [Origin_ID{leftout} ; RandID(k_folds * EachFold_Quantity + j)];
        end
        testinds = Origin_ID{leftout};
        disp(['testind is ' num2str(testinds')])
        traininds = setdiff(RandID,testinds);
    end
    num_subjects_in_fold=length(testinds);
    % Assign Subjects_data and Subjects_score data to train and test groups 
    train_data = Subjects_data;
    train_data(:,:,testinds) = [];
    train_vcts = reshape(train_data,[],size(train_data,3));
    test_data = Subjects_data(:,:,testinds); 
    train_score = Subjects_score;
    train_score(testinds) = [];
    test_score(leftout,1:num_subjects_in_fold) = Subjects_score(testinds);   
    % correlate all edges with behavior  
    if strcmp(corr_type, 'Pearson') || notDefined(corr_type)
        %correlate all edges with behavior using Pearson correlation
        [r_mat,p_mat] = corr(train_vcts',train_score);
        r_mat = reshape(r_mat,number_node,number_node);
        p_mat = reshape(p_mat,number_node,number_node);
    elseif strcmp(corr_type, 'Spearman')
        %correlate all edges with behavior using rank correlation
        [r_mat,p_mat] = corr(train_vcts',train_score,'type','Spearman');
        r_mat = reshape(r_mat,number_node,number_node);
        p_mat = reshape(p_mat,number_node,number_node); 
    elseif strcmp(corr_type, 'Partial')
        %correlate all edges with behavior using partial correlation
        train_age = all_age;%以年龄为例做偏相关，控制年龄因素的影响
        train_age(leftout) = [];
        [r_mat,p_mat] = partialcorr(train_vcts',train_score,train_age);   
        r_mat = reshape(r_mat,number_node,number_node);
        p_mat = reshape(p_mat,number_node,number_node); 
    else strcmp(corr_type, 'Robust')
        %correlate all edges with behavior using robust regression
        edge_number = size(train_vcts,1);
        r_mat = zeros(1,edge_number);
        p_mat = zeros(1,edge_number);
        for edge_i = 1 : edge_number;
            [~,stats] = robustfit(train_vcts(edge_i,:)',train_score);
            cur_t = stats.t(2);
            r_mat(edge_i) = sign(cur_t)*sqrt(cur_t^2/(number_sub-1-2+cur_t^2));
            p_mat(edge_i) = 2*tadf(cur_t,number_sub-1-2);%two talied
        end    
        r_mat = reshape(r_mat,number_node,number_node);
        p_mat = reshape(p_mat,number_node,number_node);
    end
    if  notDefined('Mask_Method'); 
        % set threshold and define masks
        pos_mask = zeros(number_node,number_node);
        neg_mask = zeros(number_node,number_node);
        pos_edges = find(r_mat > 0 & p_mat < pthresh); 
        neg_edges = find(r_mat < 0 & p_mat < pthresh);
        pos_mask(pos_edges) = 1;
        neg_mask(neg_edges) = 1;
        % pos_mask = (+(r_mat > 0 & p_mat < thresh));
        % neg_mask = (+(r_mat < 0 & p_mat < thresh)); 
    else
        %------------------------sigmoidal weighting------------------------------%
        pos_edges = find(r_mat > 0);
        neg_edges = find(r_mat < 0);
        % covert p threshold to r threshold
        T = tinv(thresh/2, number_sub-1-2);
        R = sqrt(T^2/(number_sub-1-2+T^2));
        % create a weighted mask using sigmoidal function
        % weight = 0.5, when correlation = R/3;
        % weight = 0.88, when correlation = R;
        pos_mask(pos_edges) = sigmf( r_mat(pos_edges), [3/R, R/3]);
        neg_mask(neg_edges) = sigmf( r_mat(neg_edges), [-3/R, R/3]);
    %------------------------sigmoidal weighting------------------------------%
    end
    % get sum of all edges in TRAIN subs (divide by 2 to control for the
    % fact that matrices are symmetric)
    % 获得训练集各样本所有边的总和，由于矩阵是对称的，所以除以2
    train_sumpos = zeros(length(traininds),1);
    train_sumneg = zeros(length(traininds),1);
    
    for ss = 1 : size(train_sumpos);
         % 获得训练集各样本所有pos边的总和，由于矩阵是对称的，所以除以2
        train_sumpos(ss) = sum(sum(Subjects_data(:,:,ss).*pos_mask))/2;
         % 获得训练集各样本所有neg边的总和，由于矩阵是对称的，所以除以2
        train_sumneg(ss) = sum(sum(Subjects_data(:,:,ss).*neg_mask))/2;
    end
    
    if ~Model_Method;
        % build model on TRAIN subs
        fit_pos = polyfit(train_sumpos,train_score,1);
        fit_neg = polyfit(train_sumneg,train_score,1);
        % run model on TEST subs
        % For each subject, create summary feature and use model to predict y
        for i=1:size(test_data,3)
            test_sumpos(i) = sum(sum(test_data(:,:,i).*pos_mask))/2;
            test_sumneg(i) = sum(sum(test_data(:,:,i).*neg_mask))/2;
            pred_test_score_pos(testinds(i),1) = fit_pos(1)*test_sumpos(i) + fit_pos(2);
            pred_test_score_neg(testinds(i),1) = fit_neg(1)*test_sumneg(i) + fit_neg(2);
        end
    else
        % build model on TRAIN subs
        % combining both postive and negative features
        fit_com = regress(train_score,[train_sumpos, train_sumneg, ones(num_subjects-length(testinds),1)]);
        % run model on TEST sub
        for i=1:size(test_data,3)
            test_sumpos(i) = sum(sum(test_data(:,:,i).*pos_mask))/2;
            test_sumneg(i) = sum(sum(test_data(:,:,i).*neg_mask))/2;
            pred_test_score_com(testinds(i),1) = fit_com(1)*test_sumpos(i) + fit_com(2)*test_sumneg(i) + fit_com(3);
        end
    end
end
 % compare predicted and observed scores
if ~Model_Method;
    [Prediction.r_pos, Prediction.p_pos] = corr(pred_test_score_pos,Subjects_score(RandID));
    [Prediction.r_neg, Prediction.p_neg] = corr(pred_test_score_neg,Subjects_score(RandID));
    % compare predicted and observed scores using mean squared error
    Prediction.MSE_pos        = sum((pred_test_score_pos - Subjects_score(RandID)).^2)/num_subjects;
    Prediction.MSE_neg        = sum((pred_test_score_neg - Subjects_score(RandID)).^2)/num_subjects;
else
    [Prediction.r_com, Prediction.p_com] = corr(pred_test_score_com,Subjects_score(RandID));
    Prediction.MSE_com        = sum((pred_test_score_com - Subjects_score(RandID)).^2)/num_subjects;
end

 