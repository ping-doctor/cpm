function afq = AFQ_Create(varargin)
% Create AFQ structure and set default parameters
%
%    afq = AFQ_Create(varargin)
%
% Creates anautomated fiber quantification (AFQ) structure.  The default
% fields are put in place and filled with default values.  The default
% parameters can also be changed and this will affect later stages of the
% AFQ pipeline.  The arguments to the function are in the form
% 'parameter','value'.  For example, if you call the function with
%
%    AFQ_Create('cutoff',[5 95],'sub_group',[1 1 1 0 0 0]);
% The afq.params.cutoff parameter will be set to [5 95] rather than the
% default which is [10 90] and the subject groups will be defined for the
% six subjects (3 patients and 3 controls).
% 
%
% See Also: AFQ_Set AFQ_Get
%
% Examples:
%    afq = AFQ_Create;
%    afq = AFQ_Create('runmode','test');
%
%
% (c) Jason D. Yeatman December 2011
%

%% Define the type of structure
afq.type = 'afq version 1';
%% Names of all the fiber groups
afq.fgnames = {'Left Thalamic Radiation','Right Thalamic Radiation','Left Corticospinal','Right Corticospinal', 'Left Cingulum Cingulate', 'Right Cingulum Cingulate'...
    'Left Cingulum Hippocampus','Right Cingulum Hippocampus', 'Callosum Forceps Major', 'Callosum Forceps Minor'...
    'Left IFOF','Right IFOF','Left ILF','Right ILF','Left SLF','Right SLF','Left Uncinate','Right Uncinate','Left Arcuate','Right Arcuate'};

%% Names of the ROIs associated with each fiber group
afq.roi1names = {'ATR_roi1_L', 'ATR_roi1_R','CST_roi1_L','CST_roi1_R','CGC_roi1_L.mat','CGC_roi1_R.mat'...
    'HCC_roi1_L','HCC_roi1_R','FP_L','FA_L','IFO_roi1_L','IFO_roi1_R','ILF_roi1_L','ILF_roi1_R'...
    'SLF_roi1_L','SLF_roi1_R','UNC_roi1_L','UNC_roi1_R','SLF_roi1_L','SLF_roi1_R'};
afq.roi2names = {'ATR_roi2_L', 'ATR_roi2_R','CST_roi2_L','CST_roi2_R','CGC_roi2_L.mat','CGC_roi2_R.mat'...
    'HCC_roi2_L','HCC_roi2_R','FP_R','FA_R','IFO_roi2_L','IFO_roi2_R','ILF_roi2_L','ILF_roi2_R'...
    'SLF_roi2_L','SLF_roi2_R','UNC_roi2_L','UNC_roi2_R','SLFt_roi2_L','SLFt_roi2_R'};

%% Attach a structure of values to the afq structure
vals.fa = {};
vals.md = {};
vals.rd = {};
vals.ad = {};
vals.cl = {};
afq.vals = vals;

%% Attach a cell array of subject names to the afq structure
afq.sub_names = {};

%% Attach a vector of subject groups to afq structure
afq.sub_group = [];

%% Attach a vector of subject id numbers to the afq structure
afq.sub_nums = [];

%% Attach a cell array of subject directories to the afq structure
afq.sub_dirs = {};

%% Attach a vector to define each subject's group
afq.sub_group = [];

%% Attach the tract profile structure to the afq structure

afq.TractProfiles = AFQ_CreateTractProfile;
% % Add a field to afq.TractProfiles for each tract defined above
% for ii = 1:length(fgNames)
%     afq.TractProfiles.(fgNames{ii}) = AFQ_CreateTractProfile('name',fgNames{ii});
% end

%% Attatch a field for spatial normalization
afq.xform.sn = [];
afq.xform.invDef = [];
afq.xform.ants = [];
afq.xform.antsinv = [];

%% Check which software packages are installed
afq.software.mrvista = check_mrvista;
afq.software.mrtrix = check_mrTrix;
afq.software.spm = check_spm;
afq.software.ants = check_ants;

%% Set the afq.params structure with default parameters
%  cutoff: The percentile cutoff to be used to determine what is "abnormal"
%  The default is cutoff = [10 90] meaning that subjects who fall below the
%  10th percentile or above the 90th percentile will be considered
%  "abnormal" and their data will be plotted with respect to the control
%  population
%  Example: params.cutoff = [10 90];
afq.params.cutoff = [10 90];
%  run_mode: normally leave this blank unless you just want to test this
%  script on a new dataset.  If this is the case then set run_mode = 'test'
%  in which case many fewer fibers will be tracked and segmented to allow
%  for fast testing/debugging of the script
afq.params.run_mode = [];
% cleanFibers = 1 means that once all the fiber groups have been segmented
% they will be cleaned such that any fiber that (1) any fiber that is more
% than params.maxLen standard deviations above the mean fiber length will
% be removed and and (2) any fiber that is more than params.maxDist
% standard deviations from the core of the tract will be removed.  this
% means that fibers groups will be forced to be a compact bundle
afq.params.cleanFibers = 1;
% Maximum number of iteration of the cleaning algorithm
afq.params.cleanIter = 5;
% Remove fibers that are more than maxDist standard deviations from the
% core of the tract
afq.params.maxDist = 5;
% Remove fibers that are more than maxLen standard deviations from the
% mean fiber length
afq.params.maxLen = 4;
% The number of nodes to represent each fiber
afq.params.numberOfNodes = 100;
% Should we analyze the whole length of the fiber group or just the central
% portion spanning between 2 ROIs
afq.params.clip2rois = 1;
% Set the amount of weighting that will be applied when calculating tract
% profiles. 0 meanse that each fiber contributes equally. 1 means that we
% apply gaussian weighting where each fibers contribution to the
% measurement at a node is weighted by its gaussian distance from the tract
% core. Values greater than 1 mean that the core of the tract is weighted
% more heavily and fibers futher from the core are weighted less heavily
% (faster fall off than a gaussian).  See AFQ_ComputeTractProperties
afq.params.fiberWeighting = 1;
% If params.cleanFibers==1, then this will indicate whether to perform the
% cleaning on just the clipped portion of the tract or on the full tract.
% This may be helpful for tracts like the ILF that are extremely messy with
% looping fibers etc.
afq.params.cleanClippedFibers = 0;
% The directory to save all output figures and results
afq.params.outdir = [];
% The name to give the saved afq structure
afq.params.outname = [];
% Show figures? yes or no
afq.params.showfigs = true;
% Save figures? yes or no
afq.params.savefigs = false;
% Whether or not to compute constrained spherical deconvolution using
% mrtrix
afq.params.computeCSD = 0;
% Whether or not to comput control group norms
afq.params.computenorms = 1;
% Which software package to use for normalization
afq.params.normalization = 'spm';

%% AFQ Fiber Tracking parameters
% Do fiber tracking with mrdiffusion by default. The other option is
% 'mrtrix' if it is installed and the data is HARDI
afq.params.track.algorithm = 'mrdiffusion';
% Distance between steps in the tractography algoithm
afq.params.track.stepSizeMm = 1;
% Stopping criteria FA<0.2
afq.params.track.faThresh = 0.2;
% Discard Fibers shorter than 50mm or longer than 250mm
afq.params.track.lengthThreshMm = [50 250];
% Stopping criteria angle between steps >30 degrees
afq.params.track.angleThresh = 30;
% Unknown.....
afq.params.track.wPuncture = 0.2;
% There are multip.e algorithms that can be used.  We prefer STT. See:
% Basser PJ, Pajevic S, Pierpaoli C, Duda J, Aldroubi A. 2000.
% In vivo fiber tractography using DT-MRI data.
% Magnetic Resonance in Medicine 44(4):625-32.
afq.params.track.whichAlgorithm = 1;
% Interpolation method. After each step we interpolate the tensor at that
% point. Trilinear interpolation works well.
afq.params.track.whichInterp = 1;
% This adds some randomness to each seed point. Each seed point is move
% randomly by randn*.1mm
afq.params.track.offsetJitter = 0;
% We seed in voxel in multiple locations. [0.25 and 0.75] Seeds each voxel
% at 8 equidistant locations.  For one seed in the middle of the voxel use
% afq.params.track.seedVoxelOffsets = 0.5;
afq.params.track.seedVoxelOffsets = [0.25 0.75];
% Mask from which to initialize tracking
afq.params.track.faMaskThresh = 0.30;

% TODO:
%  Write a parameter translation routine based on mrvParamFormat()
%  This will translate all of the parameters into the right format for the
%  afq parameters structure.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modify default parameters based on user input                          %
afq = afqVarargin(afq, varargin);                                         %
afq.params = afqVarargin(afq.params, varargin);     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Modify tracking parameters if the mode is test mode
if strcmp(afq.params.run_mode,'test')
    afq.params.track.seedVoxelOffsets = [0.5];
    afq.params.track.faMaskThresh = 0.35;
end
%% Attach a structure pointing to each subjects data files
for ii = 1: AFQ_get(afq,'num subs')
    temp_dtFile = ls([afq.sub_dirs{ii} '\*dt6.mat']);
    afq.files.dt6{ii} = fullfile(afq.sub_dirs{ii},temp_dtFile);
%     afq.files.dt6{ii} = fullfile(afq.sub_dirs{ii},'dt6.mat');
    if ~exist(AFQ_get(afq,'dt6 path',ii))
        error('%s file does not exist',AFQ_get(afq,'dt6 path',ii))
    end
end
afq.files.images            = struct('name',{},'path',{});
afq.files.fibers.wholebrain = cell(AFQ_get(afq,'num subs'),1);
afq.files.fibers.segmented  = cell(AFQ_get(afq,'num subs'),1);
afq.files.fibers.clean      = cell(AFQ_get(afq,'num subs'),1);

%% Add files from previous AFQ runs to afq structure

% The name of the segmented fiber group depends on whether we are clipping
% it to the ROIs or not. Or it can be passed in by the user
s = strcmp('segname',varargin) + strcmp('segName',varargin);
if sum(s)>0
    segName = varargin{find(s)+1};
elseif AFQ_get(afq,'clip2rois') == 0
    segName = 'MoriGroups_Cortex.mat';
else
    segName = 'MoriGroups.mat';
end
for ii = 1:AFQ_get(afq,'num subs')
    fibDir = fullfile(afq.sub_dirs{ii},'fibers');
    wholebrainFG = fullfile(fibDir,'WholeBrainFG.mat');
    segmentedFG = fullfile(fibDir,segName);
    cleanFG = fullfile(fibDir,[prefix(segName) '_clean_D' num2str(afq.params.maxDist) '_L'  num2str(afq.params.maxLen) '.mat']);
    if exist(wholebrainFG,'file')
        afq.files.fibers.wholebrain{ii} = wholebrainFG;
    end
    if exist(segmentedFG,'file')
        afq.files.fibers.segmented{ii} = segmentedFG;
    end
    if exist(cleanFG,'file')
        afq.files.fibers.clean{ii} = cleanFG;
    end
end
% Save the name name of the segmented fiber group
afq.files.fibers.segName = segName;
%% Allow previous analyses to be overwritten
afq.overwrite.fibers.wholebrain = zeros(AFQ_get(afq,'num subs'),1);
afq.overwrite.fibers.segmented = zeros(AFQ_get(afq,'num subs'),1);
afq.overwrite.fibers.clean = zeros(AFQ_get(afq,'num subs'),1);
afq.overwrite.vals = zeros(AFQ_get(afq,'num subs'),1);

%% If desired compute constrained spherical deconvolution with mr trix

% If mr trix is installed and CSD is computed then perform tracking on
% constrained spherical deconvolution
if afq.software.mrtrix == 1 && afq.params.computeCSD == 1
    afq.params.track.algorithm = 'mrtrix';
end

if AFQ_get(afq,'use mrtrix')
    for ii = 1:AFQ_get(afq,'num subs')
        files = AFQ_mrtrixInit(AFQ_get(afq, 'dt6 path',ii));
        afq.files.mrtrix.csd{ii} = files.csd;
        afq.files.mrtrix.wm{ii} = files.wm;
    end
end

%% Set the current subject field to subject 1
afq.currentsub = 1;

%% Add a field for meta data (eg. age, sex etc.)
afq.metadata = [];

return
