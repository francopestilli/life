function [dt6_noOR fe] = v_lifeConnectomePredictDiffusionResidual_hiromasa
%
% Generate a DWI volume with the residual signal from the fit of a Connectome.
% At every location fo where the connectome passes through.
% 
% v_lifeConnectomePredictDiffusionResidual_hiromasa
%
% Franco (C) 2012 Stanford VISTA team.

%% Initialize a connectome
baseDir    = '/biac2/wandell2/data/diffusion/pestilli/20120718_2975';
dtFile     = fullfile(baseDir,'96dirs_b2000_1point5iso_1','dt6.mat');
dwiDir     = fullfile(baseDir,'preprocessed');
dwiFile    = fullfile(dwiDir,'run01_fliprot_aligned_trilin.nii.gz');
bvecsFile  = fullfile(dwiDir,'run01_fliprot_aligned_trilin.bvecs');
bvalsFile  = fullfile(dwiDir,'run01_fliprot_aligned_trilin.bvals');
t1FileName = fullfile(baseDir,'96dirs_b2000_1point5iso_1','t1','t1.nii.gz');
fgFileName = fullfile(baseDir,'fibers','conTrack','OR-RH_FromMrTrix-20120921.pdb');
fitType    = 'voxelwise'; % Choose whether to us eth fit of life done voxel-wise or gloabally.

% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

% These are the two output files from this script that you will use with
% contrack
outDwiFile = fullfile(baseDir,sprintf('dti_noOR_%s', fitType),'raw','dwi_OR_mrtrix.nii.gz');
outBaseDir = fullfile(baseDir,sprintf('dti_noOR_%s', fitType));
mkdir(outBaseDir); % Make the directory for the new dt6 file, the one w/o OR

%% Create a Connectome
fe = feConnectomeInit(dwiFile,dtFile,fgFileName);

%% Build a DWI volume containing the residual signal from the connectome

switch fitType
  case {'voxelwise'}
    %% Fit the model one voxel at the time.
    fe = feFitModelByVoxel(fe);
    resSig = reshape(feGet(fe,'res sig full voxfit'),feGet(fe,'nbvecs'),feGet(fe,'nvoxels'));

  case {'global'}
    %% Fit the model with one weight per fiber, global weights
    fefit  = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
    fe     = feSet(fe,'fit',fefit);
    
    % Generate the residual signal
    resSig = feGet(fe,'fiber res sig with mean voxel');
  otherwise
    keyboard
end
resSig(resSig < 0) = 0; % return only positive residuals

% Load the original DWI data nifti file in which we will replace the signal
% at each location where the connectome makes a prediction
pDwi = readFileNifti(feGet(fe,'dwi file'));

% Add the Residual signal into pDwi
coords = feGet(fe,'roi coords');
%mSig = resSig;
for ic = 1:size(coords,1)
  xyz = coords(ic,:);
  %mSig(:,ic) = pDwi.data(xyz(1),xyz(2),xyz(3),feGet(fe,'bvecs indices'));
  %fprintf('Measured:%2.2f - Res%2.2f\n',mSig(5,ic),int16(resSig(5,ic)))
  pDwi.data(xyz(1),xyz(2),xyz(3),feGet(fe,'bvecs indices')) = int16(resSig(:,ic));
end

% Save the pDWi file
pDwi.fname = outDwiFile;
niftiWriteMatlab(pDwi)

%% Recompute the tensors with the new data. This operation will generate a new dt6 file.
nBootStraps = 500; % This will generate a dispersion image also (the reliability of the tesnro model fit via bootstrap)
dt6_noOR = dtiRawFitTensorMex(pDwi, bvecsFile, bvalsFile, outBaseDir,nBootStraps,[], [],[],[],1);

% Add the bvecs, bvals and t1-files to the dt6 file just created.
% Add the raw data file names (these will be full paths)
files.alignedDwRaw   = outDwiFile;
files.alignedDwBvecs = bvecsFile;
files.alignedDwBvals = bvalsFile;
save(dt6_noOR,'files','-APPEND');

return