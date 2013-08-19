function s_lifeTestTractographyModels(trackingType)
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

if notDefined('trackingType'),trackingType = 'deterministic';end

%% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

%% Start loading the data
% Basic directories
projectDir     = 'LiFE_TEST_MRTRIX_PARAMS';
basePath      = fullfile('/biac2','wandell6');
dataRootPath  = fullfile('/biac2','wandell6','data','frk','life_dti','FP20120420');
subfolders    = fullfile('150dirs_b1000_1');
baseDir       = fullfile(dataRootPath,subfolders);
dtFile        = fullfile(baseDir,'dt6.mat');
dwiFile       = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
dwiFileRepeat = fullfile(dataRootPath,'raw','0011_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
testVol       = 'V1234LO12_roi.mat';
saveDir       = fullfile(baseDir,projectDir);
t1File        = fullfile(dataRootPath,'t1','t1.nii.gz');
%trackingType  = 'deterministic'; 
switch trackingType
  case {'deterministic'}
    connectomeFile = { ...
      'mrtrix_whole_brain_lmax2_stream_ntracks_200000.pdb', ...
      'mrtrix_whole_brain_lmax2_stream_ntracks_200000_2.pdb', ...
      'mrtrix_whole_brain_lmax4_stream_ntracks_200000.pdb',...
      'mrtrix_whole_brain_lmax4_stream_ntracks_200000_2.pdb', ...
      'mrtrix_whole_brain_lmax8_stream_ntracks_200000.pdb', ...
      'mrtrix_whole_brain_lmax8_stream_ntracks_200000_2.pdb' ...
      'mrtrix_whole_brain_lmax16_stream_ntracks_200000.pdb', ...
      'mrtrix_whole_brain_lmax16_stream_ntracks_200000_2.pdb', ...
      'fg_whole_brain_FACT_RK4.pdb', ...
      'fg_whole_brain_TEND_RK4.pdb'};
  case {'probabilistic'}
    connectomeFile = { ...
      'mrtrix_whole_brain_lmax2_prob_ntracks_200000.pdb', ...
      'mrtrix_whole_brain_lmax2_prob_ntracks_200000_2.pdb', ...
      'mrtrix_whole_brain_lmax4_prob_ntracks_200000.pdb', ...
      'mrtrix_whole_brain_lmax4_prob_ntracks_200000_2.pdb', ...
      'mrtrix_whole_brain_lmax8_prob_ntracks_200000.pdb', ...
      'mrtrix_whole_brain_lmax8_prob_ntracks_200000_2.pdb' ...
      'mrtrix_whole_brain_lmax16_prob_ntracks_200000.pdb', ...
      'mrtrix_whole_brain_lmax16_prob_ntracks_200000_2.pdb'};
  otherwise
    keyboard
end

% Load the ROI defining the volume we will use to test the connectome.
volRoi = dtiReadRoi(fullfile(baseDir,projectDir,'rois',testVol));
maxVolDist = 2.5;% Max distance in mm from the ROI edges. 
    
for ii = 1:length(connectomeFile)
  fprintf('Processing: %s\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n',connectomeFile{ii})
  
  %% Set up file names
  % mapSaveNameGlobal = fullfile(saveDir,'maps',sprintf('%s_r2z_global',connectomeFile{ii}));
  % mapSaveNameVoxels = fullfile(saveDir,'maps',sprintf('%s_r2z_voxelwise',connectomeFile{ii}));
  % fgSaveName = fullfile(saveDir,'fg',sprintf('%s_clippedFG',connectomeFile{ii}));
  feSaveName = [connectomeFile{ii}(1:end-4)];
  
  %% Load a full-brain Connectome that was preprocessed with...
  disp('Loading the whole brain connectome...')
  wholeBrainPath = fullfile(dataRootPath,'fe_fibers','whole_brain_connectomes_pdb',connectomeFile{ii});
  fg = fgRead(wholeBrainPath);
  
  %% Clip the connectome to be constrained within the volumeRoi.
  tic, fprintf('Clipping Connectome fibers that leave the volume ROI...\n');
  fg = feClipFibersToVolume(fg,volRoi.coords,maxVolDist);
  fprintf('process completed in %2.3fhours\n',toc/60/60);
  
  %% Build a model of the connectome.
  tensorModelParams = [1,0];
  fe = feConnectomeInit(dwiFile,dtFile,fg,feSaveName,saveDir,dwiFileRepeat,t1File,tensorModelParams);
  clear fg
  
  %% Fit the model with global weights.
  fefit = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsig demeaned'),'sgdnn');
  fe    = feSet(fe,'fit',fefit);
  
  %% Fit the model with voxel-wise weights.
  fe = feFitModelByVoxel(fe);
  
 % %% Compute maps of R2 for each fit.
 % feSaveMapToNifti(fe,  'vox r2 zero ', mapSaveNameGlobal );
 % feSaveMapToNifti(fe,  'vox r2 zero', mapSaveNameVoxels );
  
  %% Save all files
  %fgWrite(fg, fgSaveName,'pdb');
  feConnectomeSave(fe,'_stick1')
  clear fe  
  fprintf('DONE processing: %s\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n',connectomeFile{ii})
end

% Handling paralle processing.
if ~poolwasopen, matlabpool close; end

