function fe = v_life
%
% Validates a series of calls to the LiFE functions.
% 
% v_life
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

%% Initialize the Connectome
fe         = feConnectomeInit(dwiFile,dtFile,fgFileName);

%% Estimate the weights and install them in the fe structure
fefit      = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
fe         = feSet(fe,'fit',fefit);

%% Now cross-validate the quality fo fit and install the result in the fe structure
fexval = feXvalidate(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
fe         = feSet(fe,'xvalfit',fexval);

%% Get the voxel coordinates
indexes = [2 123 1000 50];
coords  = feGet(fe,'roi coords');
coords  = coords(indexes,:); 

%% Get back the indexes from the coordinates
voxelsIndex = feGet(fe,'find voxels',coords);
assertEqual(sort(indexes(:)),find(voxelsIndex))

%% Pass coordinates
dwiVoxel = feGet(fe,'dsiinvox',coords);
dwiVoxel = feGet(fe,'dsiinvoxdemeaned',coords);   
dSig = feGet(fe,'dSig full',coords);
wiso = feGet(fe,'iso weights',coords);
res  = feGet(fe,'resfibermeanvox',coords);
res  = feGet(fe,'resfibervox',coords);
res  = feGet(fe,'residualsignalfullvoxel',coords);
res  = feGet(fe,'voxelrmse',coords);
res  = feGet(fe,'dsigfullvox',coords);
res  = feGet(fe,'varexpvox',coords);

%% Pass voxel indexes
dwiVoxel   = feGet(fe,'dsiinvox',indexes);
dwiVoxel   = feGet(fe,'dsiinvoxdemeaned',indexes);  
dSig       = feGet(fe,'dSig full',indexes);
wiso       = feGet(fe,'iso weights',indexes);
res        = feGet(fe,'resfibermeanvox',indexes);
res        = feGet(fe,'resfibervox',indexes);
res2       = feGet(fe,'residualsignalfullvoxel',indexes);
res2       = feGet(fe,'voxelrmse',indexes);
res2       = feGet(fe,'dsigfullvox',indexes);
res2       = feGet(fe,'varexpvox',indexes);


%% Fibers
fiberList = [1 100];
fg       = feGet(fe,'fiberssubset',fiberList);

%% Basic returns
fg    = feGet(fe,'fg acpc');
bvecs = feGet(fe,'bvecs');
bvals = feGet(fe,'bvals');

