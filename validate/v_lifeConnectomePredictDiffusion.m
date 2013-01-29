function fe = v_lifeConnectomePredictDiffusion
%
% Generate a synthetic volume with the signal predicted by a Connectome.
% Illustrate how to use the fit of a Connectome to generate a prediction of signal in the DWI volume.
% 
% v_lifeConnectomePredictDiffusion
%
% Franco (C) 2012 Stanford VISTA team.


%% Initialize a connectome
baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
dtFile     = fullfile(baseDir,'dti40','dt6.mat');
dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
fgFileName = fullfile(baseDir,'fibers','leftArcuateSmall.pdb');

%% Initialize the Connectome
fe     = feConnectomeInit(dwiFile,dtFile,fgFileName);

%% Estimate the weights and install them in the fe structure
fefit  = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
fe     = feSet(fe,'fit',fefit);

%% Now cross-validate the quality fo fit and install the result in the fe structure
fexval = feXvalidate(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
fe     = feSet(fe,'xvalfit',fexval);

%% Build a DWI volume containing the predicted signal from the connectome
% PRedict signal wit the connectome
pSig = feGet(fe,'psigfullvox');
dSig = feGet(fe,'dsigvox');

% Load the original DWi data nifti file in which we will replace the signal
% at each location where the connectome makes a prediction
pDwi = readFileNifti(feGet(fe,'dwi file'));

% Add the predicted signal into pDwi
coords = feGet(fe,'roi coords');
for ic = 1:size(coords,1)
  xyz = coords(ic,:);  
  %fprintf('%i DWIs: %2.3f, pSIG: %2.3f.\n',ic,(pDwi.data(xyz(1),xyz(2),xyz(3),2)),pSig(1,ic))
  pDwi.data(xyz(1),xyz(2),xyz(3),feGet(fe,'bvecsindices')) = int16(pSig(:,ic));
end

% Save the pDWi file
pDwi.fname = fullfile(baseDir,'raw','test_predicted_signal_leftArcuateSmall.nii.gz');
niftiWriteMatlab(pDwi)

% At this point a nifti file is written to disk. To use the nifti file to
% fit the same connectome used to generate it you will need to copy the
% bvecs/bvals files of the original dwi.gz.nii file into new bvecs/bvals
% files with the same file name of the newly geenrated synthetic dataset.

return