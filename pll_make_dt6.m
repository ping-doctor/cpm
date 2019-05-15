clc;
clear all;
sub_dir = 'F:\StudyData\AFQ';
sub_name = dir(sub_dir);
for i = 1 : length(sub_name)
    %if sub_name(i) is not a dir, skip
    if (isequal(sub_name(i).name, '.') == 1 || isequal(sub_name(i).name, '..') == 1)
    continue;
    end
    b0FilePath = fullfile(sub_dir, sub_name(i).name, 'native_space');
    temp_b0FileName = ls([b0FilePath '\*S0*']);
    b0FileName = fullfile(b0FilePath, temp_b0FileName);
    if ~exist('t1FileName','var')
        t1FileName = b0FileName;
        disp(t1FileName);
    end
    dt6 = dtiMakeDt6FromFsl(b0FileName, t1FileName);
end
