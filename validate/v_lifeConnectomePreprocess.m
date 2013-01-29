%% v_lifeConnectomePreprocess.m
%
% Illustrate how to preprocess a connectome (fg) to be constrained within a
% region of interest and within the cortex.
%
% We will show how to clip a fiber grou so that it only contains fibers
% that are within the volume defined by an ROI, that are not short and that
% start and end in cortex.
%
% For this example we use the volume defined by the connectome to constrain
% the fibers.
%
% See also:
%
% Franco (C) 2012 Stanford VISTA team.

%% Start loading the data
%
% At some point, we will put these fe structures into VISTADATA.
% For now, you can compute them, leave them in the feGet(fe,'savedir'), and load them
% up for testing.

% Basic directories
if ispc,  basePath = fullfile('\\red.stanford.edu','biac2-wandell6');
else      basePath = fullfile('/biac2','wandell6');
end
dataRootPath  = fullfile(basePath,'data','frk','life_dti','FP20120420');
subfolders    = fullfile('150dirs_b1000_1');
baseDir       = fullfile(dataRootPath,subfolders);
saveDir       = fullfile(baseDir,'LiFE');
dtFile        = fullfile(baseDir,'dt6.mat');
dwiFile       = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
fgFileName    = fullfile(baseDir,'fibers','whole_brain_MRTRIX','WholeBrainFG_MRTRIX.mat');
% This is the fiber group that I will use to generate the Volume ROI to
% constrain the fibers, I am working on the optic radiation so I load the
% occipital fibergroup
fgRoiName     = fullfile(saveDir,'fg','right_hemisphere_occipital_MRTRIX_FG.mat'); 
cortexFile    = fullfile(dataRootPath,'t1','t1_class.nii.gz');

% Load the fiber group to geenrat ethe Volume ROI
disp('Loading the occipital lobe fiber group...')
fg = dtiLoadFiberGroup(fgRoiName);

% Generate a volume ROI from the fiber group.
tic
fprintf('[%s] Generating an ROI for the occipital lobe...\n',mfilename)
volRoi        = dtiNewRoi('volumeRoi',rand(3,1),[]);
volRoi.coords = fefgGet(fg,'unique acpc coords');
volRoi        = dtiRoiClean(volRoi);
dtiWriteRoi(volRoi,fullfile(saveDir,'rois','volumeRoi'));
clear fg;
fprintf('done in %2.3f. ROI size %ix%i\n',toc,size(volRoi.coords))
keyboard
% Load the fiber group
disp('Loading the whole brain fiber group...')
fg = dtiLoadFiberGroup(fgFileName);

% Now restric the fibers to be withint he White-matter and the ROI
fg = feConnectomePreprocess(fg, volRoi, cortexFile, [5,6]);
fg.name = 'WholeBrainFG_MRTRIX_preprocessed';

%% Initialize the Connectome
fe     = feConnectomeInit(dwiFile,dtFile,fg);

keyboard