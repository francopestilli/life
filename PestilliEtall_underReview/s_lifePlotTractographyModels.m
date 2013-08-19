%% s_lifePlotTractographyModels.m
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
testVol       = 'V1234LO12_roi.mat';
saveDir       = fullfile(baseDir,projectDir);
trackingType  = 'deterministic'; 
switch trackingType
  case {'deterministic'}
    connectomeFile = { ...
      'mrtrix_whole_brain_lmax2_stream_ntracks_200000.mat', ...
      'mrtrix_whole_brain_lmax2_stream_ntracks_200000_2.mat', ...
      'mrtrix_whole_brain_lmax4_stream_ntracks_200000.mat',...
      'mrtrix_whole_brain_lmax4_stream_ntracks_200000_2.mat', ...
      'mrtrix_whole_brain_lmax8_stream_ntracks_200000.mat', ...
      'mrtrix_whole_brain_lmax8_stream_ntracks_200000_2.mat' ...
      'mrtrix_whole_brain_lmax16_stream_ntracks_200000.mat', ...
      'mrtrix_whole_brain_lmax16_stream_ntracks_200000_2.mat'};
  case {'probabilistic'}
    connectomeFile = { ...
      'mrtrix_whole_brain_lmax2_prob_ntracks_200000.mat', ...
      'mrtrix_whole_brain_lmax2_prob_ntracks_200000_2.mat', ...
      'mrtrix_whole_brain_lmax4_prob_ntracks_200000.mat',...
      'mrtrix_whole_brain_lmax4_prob_ntracks_200000_2.mat', ...
      'mrtrix_whole_brain_lmax8_prob_ntracks_200000.mat', ...
      'mrtrix_whole_brain_lmax8_prob_ntracks_200000_2.mat' ...
      'mrtrix_whole_brain_lmax16_prob_ntracks_200000.mat', ...
      'mrtrix_whole_brain_lmax16_prob_ntracks_200000_2.mat'};
  otherwise
    keyboard
end

for ii = 1:length(connectomeFile)
  fprintf('Processing: %s\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n',connectomeFile{ii})
  
  %% Load the connectome solution
  disp('Loading the LiFE model...')
  load(fullfile(saveDir,connectomeFile{ii}));
  
  %% Plot the figergroup
  fePlot(fe,'fg')
  keyboard
  fprintf('DONE processing: %s\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n',connectomeFile{ii})
end

% Handling paralle processing.
if ~poolwasopen, matlabpool close; end

return