function s_ms_afq_test_fascicle_load_fiberdentisty_to_mesh_afq
% Example of how to make a fider density map fr om a faccile, and render it
% on a a surface.

addpath(genpath('~/git/AFQ'));

proportionFibers = 1;

% PARAMETERS
diffusionModelParams = [1,0];       % The parameters of the tensor model AD, RD
maxVolDist           = 1;           % Max distance in mm from the ROI edges.
sdCutoff             = [4 3.32 3.32];     % One per conenctome, generally smaller for smaller lmax values
clobber              = [0 0 0 0 0]; % Owerwrite all the files.

plotFG  =0; % Shows he fiber group for the arcuate.
% DIRECTORY TO LOAD FILES FROM:
% DWI data
projectDir  = '/azure/scr1/frk/150dirs_b1000_b2000_b4000';
connectSubfolders    = {'life_mrtrix_rep1'};
% Make some plots
colors = {[.85 .5 .65],[.6 .85 .4],[.45 .6 .95]};
if notDefined('fig_saveDir'), 
    fig_saveDir = fullfile('/home/frk/Dropbox','fascicles_movies_and_figures');
    if ~exist(fig_saveDir,'dir');mkdir(fig_saveDir);end
end

% DIRECTORIES and FILES TO SAVE2:
saveDir              =   fullfile(projectDir,'results');
fe_structuresDir     =  'fe_arcuate_cst_test';
fas_structuresDir    =  'fas_arcuate_cst_test';
roi_saveDir          =  'roi_arcuate';
arcuateRoiFileName   =  'arcuate_roi_lmax12_prob';  

% Prepare the path for the nifti to the0 segemntation file
cortex = '/biac2/wandell2/data/anatomy/pestilli/t1_class_twovalued.nii.gz';
thresh = 0.0000125; % Threshold for the overlay image
crange = [.0000125 .25]; % Color range of the overlay image
smoothingKernel      = [3 3 3];
lmax = {'Tensor','ProbLmax6','ProbLmax12'};

% Select the dwi file
bval  = [2000];
for bb = 1:length(bval)
  thisbval = bval(bb);
if thisbval == 1000
  dataRootPath  = fullfile('/biac2/wandell2/data/diffusion/pestilli/20120420_2290');
  subfolders    = fullfile('150dirs_b1000_1');
  baseDir       = fullfile(dataRootPath,subfolders);
  dtFile        = fullfile(baseDir,'dt6.mat');
  dwiFile       = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
  dwiFileRepeat = fullfile(dataRootPath,'raw','0011_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
  connectomeFile= { '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax6_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_prob-500000.pdb'};
  t1File      = fullfile(dataRootPath,'t1','t1.nii.gz');

elseif thisbval == 2000
  dataRootPath  = fullfile('/biac2/wandell2/data/diffusion/pestilli/20110922_1125');
  subfolders    = fullfile('150dirs_b2000_1');
  baseDir       = fullfile(dataRootPath,subfolders);
  dtFile        = fullfile(baseDir,'dt6.mat');
  dwiFile       = fullfile(dataRootPath,'raw','0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin.nii.gz');
  dwiFileRepeat = fullfile(dataRootPath,'raw','0007_01_DTI_2mm_150dir_2x_b2000_aligned_trilin.nii.gz');
  connectomeFile= {'0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_dwi_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_tensor-500000.pdb', ...
                   '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax6_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb', ...
                   '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax12_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb'};
  t1File      = fullfile(dataRootPath,'t1','t1.nii.gz');
  
elseif thisbval == 4000
  dataRootPath  = fullfile('/biac2/wandell2/data/diffusion/pestilli/20120420_2290');
  subfolders    = fullfile('150dirs_b4000_1');
  baseDir       = fullfile(dataRootPath,subfolders);
  dtFile        = fullfile(baseDir,'dt6.mat');
  dwiFile       = fullfile(dataRootPath,'raw','0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin.nii.gz');
  dwiFileRepeat = fullfile(dataRootPath,'raw','0007_01_DWI_2mm150dir_2x_b4000_aligned_trilin.nii.gz');
  connectomeFile= {'0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_csd_lmax12_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_brainmask_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_wm_prob-500000.pdb'};
  t1File      = fullfile(dataRootPath,'t1','t1.nii.gz');
else
  keyboard
end

for irep = 1:length(connectSubfolders)
  for ii = 1:length(connectomeFile)
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
 
    % Check if the folders need to be created
    checkFolders({feSaveDir,fasSaveDir,roiSaveDir});
    
    fprintf('[%s] Analyzing the following files:\n',mfilename);
    fprintf('FAS        [%s] \n',fasFullPath);
    fprintf('CleanedFAS [%s] \n',fasCleanedFullPath);
    fprintf('ROI        [%s] \n',roiFullPath);
    fprintf('Arcuate FE [%s] \n',feSaveNameArcuate);
    fprintf('Arcute+CST [%s] \n',feSaveNameArcCST);

    
  %% Check if the these fascicles were computed already
    if ~(exist([fasFullPath,'.mat'],'file') == 2) || clobber(1)
      fprintf('[%s] Creating Fascicles \n',mfilename)
      
      fprintf('[%s] Working on: \n\n%s\n\n',wholeBrainConnectome)
      % Segment the fiber group using AFQ.
      [fascicles,classification,fg] = feAfqSegment(dtFile, wholeBrainConnectome);
      
      % Save the fascicles and the cleaned whole-brain connectome
      save([fasFullPath,'.mat'],'fascicles','classification','fg','-v7.3');
      fprintf('[%s] DONE creating fascicles for\n\n%s\n\n',mfilename, wholeBrainConnectome)
    else
      fprintf('[%s] FOUND Fascicles FILE NOT PROCESSING, Loading it\n',mfilename)
      load([fasFullPath,'.mat'],'fascicles','classification');
    end
    
    %% Check if we have fascicles that were precomputed and cleaned
    if ~(exist([fasCleanedFullPath,'.mat'],'file') == 2) || clobber(2)
      fprintf('[%s] Cleaning fascicles (SD  %2.2f) \n',mfilename,sdCutoff(ii))
      
      % These are the options for cleaning the outliers from the fascicles.
      % We pass a larger Z-score for the lmax=2, so that the fascicle is
      % larger and has more likelyhood to overlap the CST
      opts.stdCutOff   = sdCutoff(ii);   % Standard deviation fo the 3D gaussian distribution used to
      % represent the fascicle when removing outliers. 3.5 z-scores
      opts.maxLen      = 20;  % Max lenght of fibers to be accepted in cm
      opts.maxNumNodes = 100; % This is used only during the computations does not actually change the nodes
      
      % Remove the outliers from the fascicles
      [fascicles, classificationCleaned] = feAfqRemoveFascicleOutliers(fascicles, classification,opts);
      
      % Save the fascicles and the cleaned whole-brain connectome
      save([fasCleanedFullPath,'.mat'],'fascicles','classificationCleaned','-v7.3');
      fprintf('[%s] DONE cleaning fascicles \n',mfilename)
    else
      fprintf('[%s] FOUND Cleaned Fascicles NOT PROCESSING, Loading it\n',mfilename)
      load([fasCleanedFullPath,'.mat'],'fascicles','classificationCleaned');
    end
    
    % oad the nifti T1 anatomy:
    anat = niftiRead(t1File);
         
    figName = sprintf('left_fascicles_ARC_sagittal_MESH_%s',lmax{ii});
    fg = fascicles(19);
    fgMakeFiberDensityNifti(fg, figName, t1File, smoothingKernel);

    % Render the cortical surface colored by the arcuate endpoint density
    AFQ_RenderCorticalSurface(cortex, 'overlay' , figName, 'crange', crange, 'thresh', thresh,'cmap','autumn');

    % Work on the figure
    set(gcf,'color','k')
    axis off
    view(-90,0);
    saveFig(gcf,fullfile(fig_saveDir,figName))

    % Arcuate and corticospinal tract 
    figName = sprintf('left_fascicles_SLF_sagittal_MESH_%s',lmax{ii});
    fg = fascicles(16);
    fgMakeFiberDensityNifti(fg, figName, t1File, smoothingKernel);
    
    % Render the cortical surface colored by the arcuate endpoint density
    AFQ_RenderCorticalSurface(cortex, 'overlay' , figName, 'crange', crange, 'thresh', thresh,'cmap','autumn');
   
    set(gcf,'color','k')
    axis off
    view(90,0); 
    lightH = camlight('right');drawnow;
    saveFig(gcf,fullfile(fig_saveDir,figName))
    close all

  end
end % Repeated tracking
end % b-value

end


%-----------------------%
function fold = checkFolders(foldersToCheck)
% Make sure the folders exist otherwise create them:
for ff = 1:length(foldersToCheck)
  [fold,fil,ext] = fileparts(foldersToCheck{ff});
  if ~( exist(fold,'dir') == 7)
    mkdir(fold)
  end
end
end


%------------------------%
function saveFig(h,figName)
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
if ~exist(fileparts(figName),'dir'),mkdir(fileparts(figName));end

eval( sprintf('print(%s,  ''-djpeg'',''-r500'' , ''-noui'', ''%s'');', num2str(h),figName));

end


% -----------------------% 
function fdNii = fgMakeFiberDensityNifti(fg, mapName, anatomyFile, smoothingKernel)
% 
% This function will extract the fiber group from an fe structure 
% and write a nifti image of fiber density in each voxel.
%
%   fdNii = feMakeFiberDensityNifti(fe)
%
% INPUTS:
%   fe              - An fe structure.
%   smoothingKernel - 3D smoothing kernel to apply to the fiber endpoint
%                     image. If set to 0 no smoothing will be applied.
%   mapName         - Full path and file name to save output image
%   best            - Allows to show all the fibers (1) or only the ones that
%                     have non-zero weight (0). Default is to show all the
%                     fibers (0).
%
% OUTPUT:
%   fdNifti -  nifti file containing the fiber density map
%
% Franco (c) Stanford Vista Team, 2013 

% The 3D smoothing kernel to apply to the nifti image so that the density
% will look more uniform
if ~exist('smoothingKernel','var') || isempty(smoothingKernel),
  smoothingKernel = [3 3 3];
end

% Build a file name if it was not passed in.
if ~exist('mapName','var') || isempty(mapName),
  mapName = [fg.name,'_fiberDensity'];
end
 
% Load the high-res anatomical file and to build a nifti image that is 
% coregistered with the segmentation and surface files.
fdNii       = niftiRead(anatomyFile);
fdNii.fname = mapName;

% Create an image of fiber density where each voxel counts the number of
% fiber endpoints
fdNii.data = dtiComputeFiberDensityNoGUI(fg, ...
                                 fdNii.qto_xyz,    ...
                                 size(fdNii.data), ...
                                 0,1,1);
clear fg

if all(smoothingKernel) > 0
  % Smooth the image with a guassian kernel
  fdNii.data = smooth3(fdNii.data,'gaussian',smoothingKernel);
end

% Write the nifti image
niftiWrite(fdNii);

end
