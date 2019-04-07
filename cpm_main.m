% Performs Connectome-Based Predictive Modeling (CPM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%         no_iterations:
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
%   Example:
%        [yhat,perf]=cpm_main(data,gF,'pthresh',0.05,'kfolds',2);
%
%   References:
%        If you use this script, please cite:
%        Shen, X., Finn, E. S., Scheinost, D., Rosenberg, M. D., Chun, M. M.,
%          Papademetris, X., & Constable, R. T. (2017). Using connectome-based
%          predictive modeling to predict individual behavior from brain connectivity.
%          Nature Protocols, 12(3), 506.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse input
sub_name = cellstr(ls('F:\StudyData\ocd\*.txt'));
for i = 1 : length(ls('F:\StudyData\ocd\*.txt'));
    Subjects_data(:,:,i) = textread(['F:\StudyData\ocd\' sub_name{i}]);
end;
Subjects_score = xlsread('F:\教程\CBM方法机器学习\behav_value.xlsx');
corr_type = 'Pearson';
%threshold for feature selection
pthresh = 0.05;
k_folds = 28;
no_iterations = 100;
Model_Method = 0
%-----------------------------------------------------------------------------------------------------------
num_subjects = size(Subjects_data,3); % get number of subjects 
number_node = size(Subjects_data,1);% get number of  nodes in the chosen brain atlas
%% Check for errors
%%[x,y]=cpm_check_errors(x,y,kfolds);
%% Train & test Connectome-Based Predictive Model
true_prediction = cpm_cv(Subjects_data,Subjects_score,pthresh,k_folds,corr_type,Model_Method)
%% Assess performance
if no_iterations
    if ~Model_Method;
        prediction_r = zeros(no_iterations,2);
        prediction_MSE = zeros(no_iterations,2);
        prediction_r(1,1) = true_prediction.r_pos;
        prediction_r(1,2) = true_prediction.r_neg;
        prediction_MSE(1,1) = true_prediction.MSE_pos;
        prediction_MSE(1,2) = true_prediction.MSE_neg;
        %creat estimate distribution of the test statistic
        %via random shuffles of data lables
        for it = 2:no_iterations
            fprintf('\n Performing iteration %d out of %d',it,no_iterations);
            new_Subjects_score   = Subjects_score(randperm(num_subjects));
            prediction           = cpm_cv(Subjects_data,new_Subjects_score,pthresh,k_folds,corr_type,Model_Method);
            prediction_r(it,1)   = prediction.r_pos;
            prediction_r(it,2)   = prediction.r_neg;
            prediction_MSE(it,1) = prediction.MSE_pos;
            prediction_MSE(it,2) = prediction.MSE_neg;         
        end
        % Corr, the bigger, the better
        pval_r_pos       = length(find(prediction_r(:,1) >= true_prediction.r_pos)) / no_iterations;
        pval_r_neg       = length(find(prediction_r(:,2) >= true_prediction.r_neg)) / no_iterations;
        % MSE, the smaller, the better
        pval_MSE_pos     = length(find(prediction_MSE(:,1) <= true_prediction.MSE_pos)) / no_iterations;
        pval_MSE_neg     = length(find(prediction_MSE(:,2) <= true_prediction.MSE_neg)) / no_iterations;
    else
        prediction_r(1,1) = true_prediction.r_com;
        prediction_MSE(1,1) = true_prediction.MSE_com;
        for it = 2:no_iterations
            fprintf('\n Performing iteration %d out of %d',it,no_iterations);
            new_Subjects_score = Subjects_score(randperm(num_subjects));
            prediction = cpm_cv(Subjects_data,new_Subjects_score,pthresh,k_folds,corr_type,Model_Method);
            prediction_r(it,1)   = prediction.r_com;
            prediction_MSE(it,1) = prediction.MSE_com;
        end
        % Corr, the bigger, the better
        pval_r_com      = length(find(prediction_r(:,1) >= true_prediction.r_com)) / no_iterations;
        % MSE, the smaller, the better
        pval_MSE_com    = length(find(prediction_MSE(:,1) <= true_prediction.MSE_com)) / no_iterations;
    end 
end
clear;
clc;
fprintf('\nDone.\n')