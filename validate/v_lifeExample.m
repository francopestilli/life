function fe = v_lifeExample
%
% Illustrate how to open up data and run a small linear fascicle evaluation
% (LIFE).
% 
% v_lifeExample
%
% We should start making fePlot(), a gateway routine that takes an fe
% structure and some other arguments and generates visualizations of
% interesting quantities.
%
%  [uData, g] = fePlot(fe,plotType,varargin);
%
%
% Franco (C) 2012 Stanford VISTA team.


%% Initialize a connectome
baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
dtFile     = fullfile(baseDir,'dti40','dt6.mat');
dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
fgFileName = fullfile(baseDir,'fibers','leftArcuateSmall.pdb');
t1FileName = fullfile(baseDir,'t1','t1.nii.gz');
repFile    = fullfile(baseDir,'raw','dwi.nii.gz');

%% Initialize the Connectome
fe     = feConnectomeInit(dwiFile,dtFile,fgFileName,mfilename,[],repFile,t1FileName);

%% Estimate the weights and install them in the fe structure
fefit  = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
fe     = feSet(fe,'fit',fefit);

%% Now cross-validate the quality fo fit and install the result in the fe structure
fexval = feXvalidate(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
fe     = feSet(fe,'xvalfit',fexval);

%% Fit the model voxel-wise, one weight per-fiber-per-voxel
fe  = feFitModelByVoxel(fe);

%% Now test a sereis of returns
% RMSE
rmsevw  = feGet(fe,'vox rmse voxel wise');
rmse    = feGet(fe,'vox rmse');
trmse   = feGet(fe,'total rmse');
trmsevw = feGet(fe,'total rmse voxel wise');

% R2
r2vw  = feGet(fe,'vox r2 voxel wise');
r2    = feGet(fe,'vox r2');
tr2vw = feGet(fe,'total r2 voxel wise');
tr2   = feGet(fe,'total r2');

% Using the repaed measure which for this example is the same data set)
% RMSE
rmsevw  = feGetRep(fe,'vox rmse voxel wise');
rmse    = feGetRep(fe,'vox rmse');
trmse   = feGetRep(fe,'total rmse');
trmsevw = feGetRep(fe,'total rmse voxel wise');

% R2
r2vw   = feGetRep(fe,'vox r2 voxel wise');
r2     = feGetRep(fe,'vox r2');
tr2vw  = feGetRep(fe,'total r2 voxel wise');
tr2    = feGetRep(fe,'total r2');
tratio = feGetRep(fe,'total rmse ratio voxel wise');
ratio  = feGetRep(fe,'vox rmse ratio');

keyboard

return

%%
% Render a nice image of the fibers using build3DFrame and
% buildSurfaceFromFrame
fePlot(fe,'connectome');



