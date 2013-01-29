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
% Basic directories
projectDir    = 'hV4connectivity_test';
basePath      = fullfile('/biac2','wandell2');
dataRootPath  = fullfile(basePath,'data','diffusion','pestilli','20120718_2975');
saveDirPath = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/150dirs_b2000_1';
subfolders    = fullfile('96dirs_b2000_1point5iso_1');
baseDir       = fullfile(dataRootPath,subfolders);
dtFile        = fullfile(baseDir,'dt6.mat');
dwiFile       = '/biac2/wandell2/data/diffusion/pestilli/20120718_2975/raw/0006_01_DWI_15mm_96dir_Nex2_Ax2b2000.nii.gz';
dwiFileRepeated= '/biac2/wandell2/data/diffusion/pestilli/20120718_2975/raw/0008_01_DWI_15mm_96dir_Nex2_Ax2b2000.nii.gz'; % This one is a second set of measurements. Fo rour case this may be 
anatomyFile   = '/biac2/wandell2/data/diffusion/pestilli/20120718_2975/t1/t1.nii.gz';

% Fibers to load
testFiberOrig  = 'FP_LH-VOF-elsewhere-20000-Clipped-LhV4-LVF13-Oct1_cleanedWithDilate-LV3d-V3AB-LVF-cleanedWithDilate.pdb';
testFpath      = fullfile(dataRootPath,'96dirs_b2000_1point5iso_cat','mrtrix_hiromasa',testFiberOrig);

%connectomeOrig = '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax10_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_prob-500000.pdb'; % FIX FIX FIX
wholeBrainPath = '/home/htakemur/HT_fMRI-DTIproject/DWI-Pestilli-HighRes/96dirs_b2000_1point5iso_1/mrtrix_hiromasa/dwi_withOR_cat_csd_dwi_withOR_cat_brainmask_dwi_withOR_cat_wm_prob-200000.pdb'; % FIX FIX FIX

% Files to save:
saveDir           = fullfile(saveDirPath,'LiFE_hiromasa',projectDir);
testVol           = '/biac2/wandell2/data/diffusion/pestilli/20120718_2975/96dirs_b2000_1point5iso_cat/ROIs/HTmodifiedROIs/VDWT/LH-VOF-volumeROI.nii.gz';
connectomeClipped = 'HT_Connectome_MRTRIX_Clipped';
feFileNameToSave  = 'fe_addition_of_VOF_to_occipital_connectome'; 

% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

%% Generate an ROI out of the fiebrs, we will evaluate the model in this Volume
volRoi = dtiRoiFromNifti(testVol,[],[testVol(end-8:end),'_ROI_VOL'],'mat',0,1);

%% Load fibers to test and clip them to the Volume ROI
disp('Loading TEST fibers...')
fgC  = fgRead(testFpath);

% Clip the Contrac fibers to volumeRoi. THis is for extra care.
maxVolDist  = 1; % mm
tic, fprintf('Clipping TEST fibers that leave the volume ROI...\n');
fgC = feClipFibersToVolume(fgC,volRoi.coords,maxVolDist);
fprintf('process completed in %2.3fhours\n',toc/60/60);

%% Load a full-brain Connectome that was preprocessed with and clip it to the vol ROI
disp('Loading the whole brain connectome...')
fg = fgRead(wholeBrainPath);

% Clip the connectome to be constrained within the volumeRoi.
tic, fprintf('Clipping Connectome fibers that leave the volume ROI...\n');
fg = feClipFibersToVolume(fg,volRoi.coords,maxVolDist);
fprintf('process completed in %2.3fhours\n',toc/60/60);

%% Combine the Connectome with the TEST fibers.
fgAll = fgMerge(fg,fgC);

%% Build a model of the connectome
% Recompute the fe structure from the original fiber group
% Initialize the Connectome
feAll = feConnectomeInit(dwiFile,dtFile,fgAll,feFileNameToSave,saveDir,dwiFileRepeated,anatomyFile);

% Fit the model
feAll = feSet(feAll,'fit',feFitModel(feGet(feAll,'Mfiber'),feGet(feAll,'dsig demeaned'),'sgdnn'));
feAll = feFitModelByVoxel(feAll);

%% Reduce the model to containing only the fibers from the connectome
% without the connection being tested.
% Save the indices of the connectome fibers (1) and of the fibers to be
% tested (0).
% connectome_indx = 1:numel(fg.fibers);
% testFibers_indx = max(connectome_indx)+1:max(connectome_indx)+numel(fgC.fibers);
% feNo            = feConnectomeSelectFibers(feAll,connectome_indx);
feNo = feConnectomeInit(dwiFile,dtFile,fg,[feFileNameToSave,'_only_tract'],saveDir,dwiFileRepeated,anatomyFile);

% Fit the model
feNo  = feSet(feNo,'fit',feFitModel(feGet(feNo,'Mfiber'),feGet(feNo,'dsig demeaned'),'sgdnn'));
feNo  = feFitModelByVoxel(feNo);

%% Save the FE stuctures
feConnectomeSave(feAll)
feConnectomeSave(feNo)

%% Compute maps of RMSE for each fit.
feSaveMapToNifti(feAll,  'vox rmse', fullfile(saveDir,'maps','rmse_ALL_fibers') );
feSaveMapToNifti(feNo,   'vox rmse', fullfile(saveDir,'maps','rmse_NO_fibers') );

%% Save all files
dtiWriteRoi(volRoi, fullfile(saveDir,'rois',volRoi.name),'mat')
fgWrite(fgC,  fullfile(saveDir,'fg','HT_mrtrix_clipped'),'pdb');
fgWrite(fg,   fullfile(saveDir,'fg','HT_whole_brain_clipped'),'pdb');
fgWrite(fgAll,fullfile(saveDir,'fg','HT_mrtrix_AND_connectome'),'pdb');

% Handling paralle processing.
if ~poolwasopen, matlabpool close; end

keyboard