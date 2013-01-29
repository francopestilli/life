function fe = v_lifeConnectomeFitVoxels
% Illustrate how to open up data and run a small linear fascicle evaluation
% (LIFE) by fitting the fibers with one weight per voxel instead of a
% global weight.
% 
% v_lifeConnectomeFitVoxels
%
%
% Franco (C) 2012 Stanford VISTA team.

%% Initialize a connectome
baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
dtFile     = fullfile(baseDir,'dti40','dt6.mat');
dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
fgFileName = fullfile(baseDir,'fibers','leftArcuate.pdb');

% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

%% Initialize the Connectome
fe     = feConnectomeInit(dwiFile,dtFile,fgFileName);

%% Fit the model one voxel at the time.
fe = feFitModelByVoxel(fe);

%% Fit with global weights
fefit = feFitModel(feGet(fe,'Mfiber'), feGet(fe,'dsigdemeaned'),'sgdnn');
fe = feSet(fe,'fit',fefit);

%% Predict the diffusion signal with the new weights
pSigVox = feGet(fe,'psigfvoxelwise');
pSigGlb = feGet(fe,'psigfiber');
mSig = feGet(fe,'dsig demeaned');


%% Show the difference in prediction
mrvNewGraphWin('Global vs. voxel-wise fit')
subplot(1,2,1),plot(mSig,pSigGlb,'ko'), 
ylabel(sprintf('Predicted diffusion signal\n(fiber component)'));
xlabel('Measured diffusion signal (demeaned)');
title('Global fit');axis equal; axis equal
subplot(1,2,2),plot(mSig,pSigVox,'ro')
xlabel('Measured diffusion signal (demeaned)');
title('Voxel-wise fit');axis equal; axis square


return
