% Save a fiber group into a quench file

%
% Uses AFQ to segment a series of connectomes.
%
% Then it extract one fascicle from a whole-brain tractography
% We use the left arcuate as an example.
% - Build and fit a LiFE model out of the fascicle.
% - adds a second fiber group, one that crosses, to the fiber group
% - builds/fits life again
% - Shows the improvement in fit across the ROI
%
% Franco (c) Stanford Vista Team 2012

% Note:
% modify val = dtiGetValFromFibers(dt6_OR_scalarImage, fiberGroup, xform, [valName='fa'], [interpMethod='trilin'], [uniqueVox=0])
%
% To get the rRMSE values and R2 values for each node in a fiber.
% Use modify feConnectomeDisplay by looking at AFQ_RenderFibers to color
% the fiber group given a specific map.

% PARAMETERS
diffusionModelParams = [1,0];       % The parameters of the tensor model AD, RD
maxVolDist           = 1;           % Max distance in mm from the ROI edges.
sdCutoff             = [3.32 4.7];     % One per conenctome, generally smaller for smaller lmax values
clobber              = [1 1 1 1 1]; % Owerwrite all the files.

% DIRECTORY TO LOAD FILES FROM:
% DWI data
projectDir  = '/azure/scr1/frk/150dirs_b1000_b2000_b4000';
connectSubfolders    = {'life_mrtrix_rep1','life_mrtrix_rep2','life_mrtrix_rep3'};

% DIRECTORIES and FILES TO SAVE2:
saveDir              =   fullfile(projectDir,'results');
fe_structuresDir     =  'fe_arcuate_cst_test';
fas_structuresDir    =  'fas_arcuate_cst_test';
roi_saveDir          =  'roi_arcuate';
arcuateRoiFileName   =  'arcuate_roi_lmax12_prob';

% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

% Select the dwi file
bval  = [1000];
bb = 1;
thisbval = bval(bb);
if thisbval == 1000
  dataRootPath  = fullfile('/biac2/wandell2/data/diffusion/pestilli/20120420_2290');
  subfolders    = fullfile('150dirs_b1000_1');
  baseDir       = fullfile(dataRootPath,subfolders);
  dtFile        = fullfile(baseDir,'dt6.mat');
  dwiFile       = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
  dwiFileRepeat = fullfile(dataRootPath,'raw','0011_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
  connectomeFile= { '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax12_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_prob-500000.pdb', ...
    '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax2_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_stream-500000.pdb'};
  t1File      = fullfile(dataRootPath,'t1','t1.nii.gz');
  
elseif thisbval == 2000
  dataRootPath  = fullfile('/biac2/wandell2/data/diffusion/pestilli/20110922_1125');
  subfolders    = fullfile('150dirs_b2000_1');
  baseDir       = fullfile(dataRootPath,subfolders);
  dtFile        = fullfile(baseDir,'dt6.mat');
  dwiFile       = fullfile(dataRootPath,'raw','0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin.nii.gz');
  dwiFileRepeat = fullfile(dataRootPath,'raw','0007_01_DTI_2mm_150dir_2x_b2000_aligned_trilin.nii.gz');
  connectomeFile= {'0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax12_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb', ...
    '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax2_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_stream-500000.pdb'};
  t1File      = fullfile(dataRootPath,'t1','t1.nii.gz');
  
elseif thisbval == 4000
  dataRootPath  = fullfile('/biac2/wandell2/data/diffusion/pestilli/20120420_2290');
  subfolders    = fullfile('150dirs_b4000_1');
  baseDir       = fullfile(dataRootPath,subfolders);
  dtFile        = fullfile(baseDir,'dt6.mat');
  dwiFile       = fullfile(dataRootPath,'raw','0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin.nii.gz');
  dwiFileRepeat = fullfile(dataRootPath,'raw','0007_01_DWI_2mm150dir_2x_b4000_aligned_trilin.nii.gz');
  connectomeFile= {'0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_csd_lmax12_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_brainmask_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_wm_prob-500000.pdb', ...
    '0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_csd_lmax2_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_brainmask_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_wm_stream-500000.pdb'};
  t1File      = fullfile(dataRootPath,'t1','t1.nii.gz');
else
  keyboard
end

irep = 1; ii = 1;
fprintf('\n\n[%s] NEW CONNECTOME ###############################################\n',mfilename);

% Name of the whole-brain connectome to load.
wholeBrainConnectome = fullfile(projectDir,connectSubfolders{irep},connectomeFile{ii});

% Set up file names
% Find an identificative name for the connectome that is shortenough:
cName = [connectomeFile{ii}(1:57),'_',connectomeFile{ii}(end-16:end-4)];

% Name and path for savign the fe structures
feSaveDir  = fullfile(saveDir,connectSubfolders{irep},fe_structuresDir);
fasSaveDir = fullfile(saveDir,connectSubfolders{irep},fas_structuresDir);
roiSaveDir = fullfile(saveDir,connectSubfolders{1},roi_saveDir);

% Build the full name of the two fascicles FE's structures
feSaveNameAll     = sprintf('%s_diffModAx%sRd%s_%s',cName,num2str(100*diffusionModelParams(1)), ...
  num2str(100*diffusionModelParams(2)));
feSaveNameArcuate = sprintf('arcuateFE_%s',feSaveNameAll);
feSaveNameArcCST  = sprintf('arcuateCstUnionFE_%s',feSaveNameAll);

% Name and path for saving the fascicles
fasSaveName        = sprintf('fas_%s',feSaveNameAll);
fasFullPath        = fullfile(fasSaveDir,fasSaveName);
fasCleanedSaveName = sprintf('%s_cleaned_sd%i',fasSaveName,100*sdCutoff(ii));
fasCleanedFullPath = fullfile(fasSaveDir,fasCleanedSaveName);

% Build a name for the roi
%arcuateRoiFileName  = sprintf('arcuateROI_%s_sd%2.0f.mat',arcuateRoiFileName,100*sdCutoff(1));
roiFullPath   = fullfile(roiSaveDir,arcuateRoiFileName);

fprintf('[%s] Analyzing the following files:\n',mfilename);
fprintf('FAS        [%s] \n',fasFullPath);
fprintf('CleanedFAS [%s] \n',fasCleanedFullPath);
fprintf('ROI        [%s] \n',roiFullPath);
fprintf('Arcuate FE [%s] \n',feSaveNameArcuate);
fprintf('Arcute+CST [%s] \n',feSaveNameArcCST);

fprintf('[%s] FOUND Cleaned Fascicles NOT PROCESSING, Loading it\n',mfilename)
load([fasCleanedFullPath,'.mat'],'fascicles','classificationCleaned');

%% LOAD DIFFERENT ARCUATE FASCICOLI To Show the differences of the fascicles with deifferent algorithms
for ii = 1:20
  quenchWriteStateFromFiberGroup(fascicles,t1File,'/Dropbox/fascicles_0%i')
end
