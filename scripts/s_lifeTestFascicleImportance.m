%% s_lifeTestFascicleImportance.m
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
% Basic directories
projectDir     = 'HT_V123_V4';
myLifeProjects = 'LiFE_hiromasa';
basePath      = fullfile('biac2','wandell6');
dataRootPath  = fullfile('biac2','wandell6','data','frk','life_dti','FP20120420');
subfolders    = fullfile('150dirs_b1000_1');
baseDir       = fullfile(dataRootPath,subfolders);
dtFile        = fullfile(baseDir,'dt6.mat');
dwiFile       = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');

% Files to load:
testFiberOrig  = 'RV123-HT-10deg.nii_1_cleaned_RV4-LVF13-Aug13_nonZero_MaskROI_20120917T151607.pdb';
testFpath      = fullfile(dataRootPath,'fibers','conTrack',testFiberOrig);

connectomeOrig = 'WholeBrainFG_MRTRIX.mat';
wholeBrainPath = fullfile(baseDir,'fibers','whole_brain_MRTRIX',connectomeOrig);

% Files to save:
saveDir           = fullfile(baseDir,myLifeProjects,projectDir);
testVol           = 'HT_RV123_RV4';
testFibersClipped = 'HT_V123_V4_Clipped';
connectomeClipped = 'HT_Connectome_Contrack_Clipped';

% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

%% Load ConTrack fibers
disp('Loading TEST fibers...')
fgC = fgRead(testFpath);

% Generate an ROI out of the fiebrs, we will evaluate the model in this
% Volume
randColor = rand(1,3);
volRoi    = dtiNewRoi(testVol,randColor,fefgGet(fgC,'unique image coords'));
volRoi    = dtiRoiClean(volRoi,3,['fillholes', 'dilate', 'removesat']);

%% Clip the TEST fibers to volumeRoi. This is for extra care.
maxVolDist  = 1; % mm
tic, fprintf('Clipping TEST fibers that leave the volume ROI...\n');
fgC = feClipFibersToVolume(fgC,volRoi.coords,maxVolDist);
fprintf('process completed in %2.3fhours\n',toc/60/60);

%% Load a full-brain Connectome that was preprocessed with
% feConnectomePreprocess.m
disp('Loading the whole brain connectome...')
fg = fgRead(wholeBrainPath);

%% Clip the connectome to be constrained within the volumeRoi.
tic, fprintf('Clipping Connectome fibers that leave the volume ROI...\n');
fg = feClipFibersToVolume(fg,volRoi.coords,maxVolDist);
fprintf('process completed in %2.3fhours\n',toc/60/60);

%% Combine the Connectome with the TEST fibers.
fgAll = fgMerge(fg,fgC);

%% Build a model of the connectome
feAll = feConnectomeInit(dwiFile,dtFile,fgAll);

% Fit the model
fefit = feFitModel(feGet(feAll,'Mfiber'),feGet(feAll,'dsig demeaned'),'sgdnn');
feAll = feSet(feAll,'fit',fefit);

%% Reduce the model to containing only the fibers from the connectome
% without the connection being tested.
% Save the indices of the connectome fibers (1) and of the fibers to be
% tested (0).
connectome_indx = 1:numel(fg.fibers);
testFibers_indx = max(connectome_indx)+1:max(connectome_indx)+numel(fgC.fibers);
feNo            = feConnectomeSelectFibers(feAll,connectome_indx);

% Fit the model
fefit = feFitModel(feGet(feNo,'Mfiber'),feGet(feNo,'dsig demeaned'),'sgdnn');
feNo  = feSet(feNo,'fit',fefit);

%% Reduce the model to containing only the TEST fibers being tested.
feTestF = feConnectomeSelectFibers(feTestF,testFibers_indx);

% Fit the model
fefit   = feFitModel(feGet(feTestF,'Mfiber'),feGet(feTestF,'dsig demeaned'),'sgdnn');
feTestF = feSet(feTestF,'fit',fefit);

%% Compute maps of RMSE for each fit.
feSaveMapToNifti(feAll,  'vox rmse', fullfile(saveDir,'maps','rmse_ALL_fibers') );
feSaveMapToNifti(feNo,   'vox rmse', fullfile(saveDir,'maps','rmse_NO_fibers') );
feSaveMapToNifti(feTestF,'vox rmse', fullfile(saveDir,'maps','rmse_test_fibers') );

%% Compute maps of R2 for each fit.
feSaveMapToNifti(feAll,  'vox r2 zero', fullfile(saveDir,'maps','r2z_ALL_fibers') );
feSaveMapToNifti(feNo,   'vox r2 zero', fullfile(saveDir,'maps','r2z_NO_fibers') );
feSaveMapToNifti(feTestF,'vox r2 zero', fullfile(saveDir,'maps','r2z_test_fibers') );

%% Save all files
dtiWriteRoi(volRoi, fullfile(saveDir,'rois',volRoi.name),'mat')
fgWrite(fgC,  fullfile(saveDir,'fg','HT_contract_clipped'),'pdb');
fgWrite(fg,   fullfile(saveDir,'fg','HT_whole_brain_clipped'),'pdb');
fgWrite(fgAll,fullfile(saveDir,'fg','HT_contrack_AND_connectome'),'pdb');
feConnectomeSave(feAll)
feConnectomeSave(feNo)
feConnectomeSave(feTestF)

% Handling paralle processing.
if ~poolwasopen, matlabpool close; end

