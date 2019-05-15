clc;
clear all;
%  Example:
% 
%   % Get the path to the AFQ directories
%   [AFQbase AFQdata] = AFQ_directories;
%   % Create a cell array where each cell is the path to a data directory
%   sub_dirs = {[AFQdata '/patient_01/dti30'], [AFQdata '/patient_02/dti30']...
%   [AFQdata '/patient_03/dti30'], [AFQdata '/control_01/dti30']...
%   [AFQdata '/control_02/dti30'], [AFQdata '/control_03/dti30']};
%   % Create a vector of 0s and 1s defining who is a patient and a control
%   sub_group = [1, 1, 1, 0, 0, 0]; 
sub_dir = 'F:\StudyData\AFQ';
sub_name = dir(sub_dir);
sub_dirs= {};
sub_group = [1, 1, 1, 0, 0, 0];
for i = 1 : length(sub_name)
    %if sub_name(i) is not a dir, skip
    if (isequal(sub_name(i).name, '.') == 1 || isequal(sub_name(i).name, '..') == 1 || sub_name(i).isdir == 0);
        continue;     
    end 
    sub_dirs(i-2) = {[sub_dir '\' sub_name(i).name '\native_space']};
end
outdir = sub_dir;
outname = fullfile(sub_dir, ['afq_' datestr(now, 'yyyy_mm_dd_HHMM')]);
afq = pll_AFQ_Create('run_mode','test', 'sub_dirs', sub_dirs, 'sub_group',sub_group, 'showfigs',false);
% Run AFQ to generate the fiber tracts
[afq patient_data control_data norms abn abnTracts] = pll_AFQ_run(sub_dirs, sub_group, afq);
%[afq, patient_data, control_data, norms, abn, abnTracts] = AFQ_run(sub_dirs,sub_group, afq);
save(outname, 'afq');
