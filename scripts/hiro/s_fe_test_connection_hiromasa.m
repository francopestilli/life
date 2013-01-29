%% s_fe_test_connection_hiromasa
%
% This file is a modification of s_fe_test_connection.m
% THe purpose is to give Hiromasa a start point for testing connections
% between different ROIs.
% 
% Hiromasa is interested in testing the connection of hV4 (upper and lower
% visual field) to V1, V2, V3, v3ab, LO1, LO2, TO1, and TO2.
%
% This script shows how to test the connection between a pair of his ROIs
%
% The logic here is:
% (1) Load a fiber group created with MRTRIX of the whole occipital lobe.
% (2) Load a fiber group created with Contrack by tracking between hV4 and
%     a cluster of ROIS (V1, V2, V3, V3ab, LO1, LO2, TO1, and TO2).
% (3) Create our Connectome: Combine the two fiber groups in (1) and (2).
% (4) Build a LiFE model (fe structure)
% (5) Cull the Connectome using feConenctomeCull.m
% (6) Save the Connectome model in a file that we will reuse several times.
% (7) Perform an 'AND' operation with the fibers and two of the ROIs of 
%     interest, for example, V1, hV4 Lower.  
%     This will create a new connectome (fe structure) let's call it
%     feNOT_v1_to_v4Lower
% (8) Compute a map of RMSE (root-mean-squared-error) with fe and
%     feNOT_v1_to_v4Lower to show how the fit changes and where it changes when
%     the connection is removed.
%
% Illustrate how to compare the quality of predictions between a connectome
% that has a full complement of fascicles, and one that has a
% specific subset of these fascicles removed.
%
% See also:
%
% Franco (c) 2012 Stanford VISTA team.

%% Start loading the data
%
% At some point, we will put these fe structures into VISTADATA.
% For now, you can compute them, leave them in the saveDir, and load them
% up for testing.

recompute_fe = 1; % if 1, it will load a fibergroup and recompute the fe structure.

% Basic directories
if ispc,  basePath = fullfile('\\red.stanford.edu','biac2-wandell6');
else      basePath = fullfile('/biac2','wandell6');
end
dataRootPath = fullfile(basePath,'data','frk','life_dti','FP20120420');
subfolders   = fullfile('150dirs_b1000_1');
baseDir      = fullfile(dataRootPath,subfolders);
dtFile       = fullfile(baseDir,'dt6.mat');
dwiFile      = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');

% Where to save the results, and the file name for the connectome:
saveDir      = fullfile(baseDir,'LiFE_hiromasa');
feFileName = 'hiromasa_hV4_project_connectome.mat';

% This fiber group was build in mrDiffusion by combining these two fiber groups:
% occipital_FG      = fullfile(saveDir,'right_occipital_MRTRIX_FG.mat')
% contract_V1_2_TO2 = PATH/TO/FILE
fgFileName           = fullfile(saveDir,'right_occipital_MRTRIX_FG_AND_contrackV1_2_TO1.mat');

% These are the ROIs to use for the test:
roisFileNames = {fullfile(saveDir,'hV4.mat'),fullfile(saveDir,'LO1.mat')};

%% Load a precomputed Connectome or Build and fit one.
if recompute_fe
  % Recompute the fe structure from the original fiber group
  % Initialize the Connectome
  fe = feConnectomeInit(dwiFile,dtFile,fgFileName);
  % Se the connectome name.
  fe = feSet(fe,'name',feFileName);
  
  % Cull the conenctome, leaving only the fibers that account for most of
  % the variance.
  [fe fefit fexval] = feCullConnectome(fe,maxNumInter, lowerPrct, prctR2redux);
  fe     = feSet(fe,'fit',fefit);
  fe     = feSet(fe,'xval fit',fexval); 
  
  % This will save inside saveDir
  feConnectomeSave(fe);
  
else
  % Load the fe structure from a precomputed file.
  disp('Loading a precomputed fe structure...')
  load(fullfile(saveDir,feFileName))
end


%% Test hypothesis: that there are important fibers connecting two ROIs
% Load a hV4 ROI
% In acpc space
hV4 = dtiReadRoi(roisFileNames{1});

% xform the ROI coordinates to IMG space
hV4.coords =  mrAnatXformCoords(feGet(fe,'xform acpc2img'),hV4.coords);

% Load a LO1 ROI
% In acpc space
LO1 = dtiReadRoi(roisFileNames{2});

% xform the ROI coordinates to IMG space
LO1.coords =  mrAnatXformCoords(feGet(fe,'xform acpc2img'),LO1.coords);

% Keep fibers touching V1
[~,~, fibersTouchingROIs] = dtiIntersectFibersWithRoi([], {'and'}, hV4, LO1, feGet(fe,'fg img'));

% Select the fibers to keep in the connectome.
% We keep the fibers that DO NOT touch the two ROIs
fibersToKeep = find(~fibersTouchingROIs);

% Reduce the connectome. We remove all the fibers touching V1.
feNOT      = feConnectomeSelectFibers(fe,fibersToKeep);

% Now cross-validate the quality fo fit and install the result in the fe structure
fefit = feFitModel(feGet(feNOT,'Mfiber'),feGet(feNOT,'dsigfiber'),'sgdnn');
feNOT = feSet(feNOTv1,'fit',fefit);

%% Create a statistical map of the quality of fit of the Original conenctome.
[sz, dwi] = dwiGet(dwiFile,'volume size');

% RMSE map
niftiName = fullfile(saveDir,'rmse_test');
rmse  = feValues2volume(feGet(fe,'vox rmse'),feGet(fe,'roi coords'),sz(1:3));
feWriteValues2nifti(rmse,niftiName,feGet(fe,'xform img 2 Acpc'));

% R2 map
niftiName = fullfile(saveDir,'r2_test');
r2  = feValues2volume(feGet(fe,'vox r2'),feGet(fe,'roi coords'),sz(1:3));
feWriteValues2nifti(r2,niftiName,feGet(fe,'xform img 2 Acpc'));

%% Create a statistical map of the quality of fit of the Connectome without the fibers touching the two ROIs
niftiName = fullfile(saveDir,'rmse_test_not');
rmse  = feValues2volume(feGet(feNOT,'vox rmse'),feGet(feNOT,'roi coords'),sz(1:3));
feWriteValues2nifti(rmse,niftiName,feGet(feNOT,'xform img 2 Acpc'));

niftiName = fullfile(saveDir,'r2_test_not');
r2  = feValues2volume(feGet(feNOT,'vox r2'),feGet(feNOT,'roi coords'),sz(1:3));
feWriteValues2nifti(r2,niftiName,feGet(feNOT,'xform img 2 Acpc'));

%% Visualize the maps in mrDiffusion: Ask Franco

%% End