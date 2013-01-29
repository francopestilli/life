%% s_fe_visualize
% 
% After computing the RMSE of two LiFE fits (see s_fe_test_connection.m),
% we visualize the results.
%
% You s
%
% Franco (c) VISTASOFT Team, 2012

%% Start loading the data
%
% At some point, we will put these fe structures into VISTADATA.
% For now, you can compute them, leave them in the tempdir, and load them
% up for testing.

% Basic directories
if ispc,  basePath = fullfile('\\red.stanford.edu','biac2-wandell6');
else      basePath = fullfile('/biac2','wandell6');
end
dataRootPath = fullfile(basePath,'data','frk','life_dti','FP20120420');
subfolders   = fullfile('150dirs_b1000_1');
baseDir      = fullfile(dataRootPath,subfolders);
saveDir      = fullfile(baseDir,'LiFE');
dtFile       = fullfile(baseDir,'dt6.mat');
dwiFile      = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
fgFileName   = fullfile(saveDir,'right_occipital_MRTRIX_FG.mat');

%% Load the fe structure from a precomputed file.
disp('Loading a precomputed fe structure...')
feFileName = 'test_hypothesis.mat';
load(fullfile(saveDir,feFileName));

%% Creating the parameter map for the overlay
[sz, dwi] = dwiGet(dwiFile,'volume size');
niftiName = fullfile(tempdir,'rmse_test');
rmse      = feValues2volume(feGet(fe,'vox rmse'),feGet(fe,'roi coords'),sz(1:3));
nii = niftiCreate('data',rmse,...
    'qto_xyz',feGet(fe,'xform img 2 acpc'),...
    'fname',niftiName);
niftiWrite(nii);

% The last two lines are equivalent to this.  FP likes this call better.
% feWriteValues2nifti(rmse,niftiName,feGet(fe,'xform img 2 Acpc'));

% If you want it in NIFTI-1 format, you can do this:
% tmp = niftiVista2ni(nii);
% Is there some niftiViewer?
%   showMOverlayImage would be a good starting point.
%
% You can check if you are worried.
% nii2 = niftiRead(niftiName);


%% Creating R2 parameter map for the overlay
niftiName = fullfile(tempdir,'r2_test');
r2        = feValues2volume(feGet(fe,'vox r2'),feGet(fe,'roi coords'),sz(1:3));
nii = niftiCreate('data',r2,...
    'qto_xyz',feGet(fe,'xform img 2 acpc'),...
    'fname',niftiName);
niftiWrite(nii);

%%  Open mrDiffusion to get the data and figure handle

% Startup up the big boy
[dtiF, dtiH] = mrDiffusion('on',dtFile);

% we need the xform that brings us inline with the anatomy image. We
% assume that the NIFTI/Analyze xform brings us to ac-pc space, so we want
% the difference between the two ac-pc xforms.
mmPerVox = nii.pixdim;
xform    = nii.qto_xyz;
unitStr  = 'err';
statName = 'rmse';
rmse = rmse/max(rmse(:));
dtiH = dtiAddBackgroundImage(dtiH, rmse, statName, mmPerVox, xform, [], 0, unitStr);

% Messing around for muscle data
%handles = dtiAddBackgroundImage(handles, img, f, mmPerVox, xform, [0.5 1.2], 0, unitStr);
dtiH = dtiRefreshFigure(dtiH);
guidata(dtiF, dtiH);

%%
mmPerVox = nii.pixdim;
xform    = nii.qto_xyz;
unitStr  = 'corr';
statName = 'r2';
r2 = r2/max(r2(:));
dtiH = dtiAddBackgroundImage(dtiH, r2, statName, mmPerVox, xform, [], 0, unitStr);

% Messing around for muscle data
%handles = dtiAddBackgroundImage(handles, img, f, mmPerVox, xform, [0.5 1.2], 0, unitStr);
dtiH = dtiRefreshFigure(dtiH);
guidata(dtiF, dtiH);

%% Load an ROI and get the rmse
%  This is in acpc space.
ROI = dtiReadRoi(fullfile(feRootPath,  'sphere_05.mat'));
ROI.coords = mrAnatXformCoords(feGet(fe,'xform acpc 2 img'),ROI.coords);
ROI.coords = unique(floor(ROI.coords),'rows');  % Perhaps this should go into the function?
rmse = feGet(fe,'voxel rmse',ROI.coords);
size(rmse)
mean(rmse)

mrvNewGraphWin; hist(rmse,20)

%% Now a 2nd ROI
ROI2 = dtiReadRoi(fullfile(feRootPath,  'sphere2_05.mat'));
ROI2.coords = mrAnatXformCoords(feGet(fe,'xform acpc 2 img'),ROI2.coords);
ROI2.coords = unique(floor(ROI2.coords),'rows');  % Perhaps this should go into the function?
rmse2 = feGet(fe,'voxel rmse',ROI2.coords);
size(rmse2)
mean(rmse2)

mrvNewGraphWin; hist(rmse2,20)

%% sphere_10
% rect_-11_-78__01
ROI3 = dtiReadRoi(fullfile(feRootPath,  'Poly_(1).mat'));
ROI3.coords = mrAnatXformCoords(feGet(fe,'xform acpc 2 img'),ROI3.coords);
ROI3.coords = unique(round(ROI3.coords),'rows');  % Perhaps this should go into the function?
rmse3 = feGet(fe,'voxel rmse',ROI3.coords);
size(rmse3)
mean(rmse3)

mrvNewGraphWin; hist(rmse3,20)

%% End