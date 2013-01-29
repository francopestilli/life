function fe = v_lifeConnectomePredictDiffusionResidual
%
% Generate a DWI volume with the residual signal from the fit of a Connectome.
% At every location fo where the connectome passes through.
% 
% v_lifeConnectomePredictDiffusionResidual
%
% Franco (C) 2012 Stanford VISTA team.


%% Initialize a connectome
baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
dwiDir     = fullfile(baseDir,'raw');
t1FileName = fullfile(baseDir,'t1','t1.nii.gz');
dtFile     = fullfile(baseDir,'dti40','dt6.mat');
dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
bvecsFile  = fullfile(baseDir,'raw','dwi.bvecs');
bvalsFile  = fullfile(baseDir,'raw','dwi.bvals');
fgFileName = fullfile(baseDir,'fibers','leftArcuateSmall.pdb');

% Name of the DWI file and of the directory for the dt6 file that will be generated at the end of
% this project
outBaseDir = fullfile(baseDir,'dti_tmp');
mkdir(outBaseDir); % Make the directory for the new dt6 file, the one w/o OR
outDwiFile = fullfile(outBaseDir,sprintf('%s.nii.gz',mfilename));

%% Initialize the Connectome
fe     = feConnectomeInit(dwiFile,dtFile,fgFileName);

%% Estimate the weights and install them in the fe structure
fefit  = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
fe     = feSet(fe,'fit',fefit);

%% Build a DWI volume containing the residual signal from the connectome
% Generate the residual signal
resSig = feGet(fe,'fiber res sig with mean voxel');
resSig(resSig < 0) = 0; % return only positive residuals

% Fit the model one voxel at the time.
%fe = feFitModelByVoxel(fe);
%resSig = feGet(fe,'res sig full voxfit');
%resSig = reshape(resSig2,feGet(fe,'nbvecs'),feGet(fe,'nvoxels'));

% Load the original DWI data nifti file in which we will replace the signal
% at each location where the connectome makes a prediction
pDwi = readFileNifti(feGet(fe,'dwi file'));

% Add the predicted signal into pDwi
coords = feGet(fe,'roi coords');
%mSig = resSig;
for ic = 1:size(coords,1)
  xyz = coords(ic,:);  
  %mSig(:,ic) = double(pDwi.data(xyz(1),xyz(2),xyz(3),feGet(fe,'bvecs indices')));
  %fprintf('%i DWIs: %2.3f %2.3f, RES: %2.3f %2.3f.\n',ic,mSig([1 80],ic),resSig([1 80],ic))
  pDwi.data(xyz(1),xyz(2),xyz(3),feGet(fe,'bvecs indices')) = int16(resSig(:,ic));
end

% Save the pDWi file
pDwi.fname = outDwiFile;
niftiWriteMatlab(pDwi)

%% Recompute the tensors with the new data. This operation will generate a new dt6 file.
nBootStraps = 10; % 500 is the default. This will generate a dispersion image also 
                  %        (the reliability of the tensor model fit via bootstrap)
dt6_tmp = dtiRawFitTensorMex(pDwi, bvecsFile, bvalsFile, outBaseDir,nBootStraps,[], [],[],[],1);

% Add the bvecs, bvals and t1-files to the dt6 file just created.
% Add the raw data file names (these will be full paths)
files.alignedDwRaw   = outDwiFile;
files.alignedDwBvecs = bvecsFile;
files.alignedDwBvals = bvalsFile;
save(dt6_tmp,'files','-APPEND');

return