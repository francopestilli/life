function s_ms_afq_fascicles_figures
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

addpath(genpath('~/git/AFQ'));

proportionFibers = 1;

% PARAMETERS
diffusionModelParams = [1,0];       % The parameters of the tensor model AD, RD
maxVolDist           = 1;           % Max distance in mm from the ROI edges.
sdCutoff             = [3.32 3.32 3.32 3.32];     % One per conenctome, generally smaller for smaller lmax values
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

% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

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
  connectomeFile= { ...
%     '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_dwi_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_tensor-500000.pdb'%, ...
     '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax6_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb'};%, ...
%     '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax10_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb'%, ...
%     '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax14_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb', ...
%     '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax16_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb'};%, ...
%     '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax2_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_stream-500000.pdb'
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

lmax = {'ProbLmax6'};

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
     
    % Arcuate and corticospinal tract 
    figName = sprintf('left_fascicles_ARC_sagittal_%s',lmax{ii});
    fg = fascicles(19);
    subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)*proportionFibers));
    fg.fibers = fg.fibers(subsampleindx);
    [figH1, lightH] = mbaDisplayConnectome(fg.fibers,figure);
    hold on
    sliceH = feDisplayBrainSlice(anat, [-15 0 0]);    
    view(-90,0);
    set(gca,'ylim',[-70 55],'xlim',[-75 0],'zlim',[-25 85])
    delete(lightH)
    lightH = camlight('right');drawnow;
    saveFig(gcf,fullfile(fig_saveDir,figName))
        
%    % The following is an attept to plot the arcuate with the FA map
    figName = sprintf('left_fascicles_ARC_sagittal_FA_map_%s',lmax{ii});
    fg = fascicles(19);
     % Eigenvalues are necessary for all the computations.
    % Here we precompute them then we pass them in for different values:
    sprintf('\nComputing eigenvalues for %s...\n',fg.name)
    eigenvals = fefgGet(fg,'eigenvals',dtFile);
    fa = fefgGet(fg, 'fa', eigenvals);
    
    subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)*proportionFibers));
    fg.fibers = fg.fibers(subsampleindx);
    faDisp = fa(subsampleindx);
    [figH1, lightH] = mbaDisplayConnectome(fg.fibers,figure,faDisp,'map',hot(60));
    hold on
    sliceH = feDisplayBrainSlice(anat, [-40 0 0]);    
    view(90,0);
    set(gca,'ylim',[-70 55],'xlim',[-75 0],'zlim',[-25 85])
    delete(lightH)
    lightH = camlight('right');drawnow;
    saveFig(gcf,fullfile(fig_saveDir,figName))
    
    % Arcuate and corticospinal tract 
    figName = sprintf('left_fascicles_ARC_CST_sagittal_FA_map_%s',lmax{ii});
    sprintf('\nComputing eigenvalues for %s...\n',fg.name)
    fgCST = fascicles(3);
    eigenvals = fefgGet(fgCST,'eigenvals',dtFile);
    faCST = fefgGet(fgCST, 'fa', eigenvals);
    
    fibers = {};
    fibers(1:length(fg.fibers)) = fg.fibers;
    fibers(length(fg.fibers)+1:length(fg.fibers)+length(fgCST.fibers)) = fgCST.fibers;
     
    faAll(1:length(fa)) = fa;
    faAll(length(fa)+1:length(fa)+length(faCST)) = faCST;

    [figH1, lightH] = mbaDisplayConnectome(fibers,figure,faAll,'map',hot(60));
    hold on
    sliceH = feDisplayBrainSlice(anat, [-40 0 0]);
    
    subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)*proportionFibers));
    fg.fibers = fg.fibers(subsampleindx);
    [figH2, lightH2] = mbaDisplayConnectome(fg.fibers,gcf,[.8 .6 .2],'single');
    view(90,0);
    set(gca,'ylim',[-70 55],'xlim',[-75 0],'zlim',[-25 85])
    delete(lightH)
    delete(lightH2)
    lightH = camlight('right');drawnow;
    saveFig(gcf,fullfile(fig_saveDir,figName))
    
    % Arcuate and corticospinal tract 
    figName = sprintf('left_fascicles_ARC_CST_sagittal_%s',lmax{ii});
    fg = fascicles(19);
    subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)*proportionFibers));
    fg.fibers = fg.fibers(subsampleindx);
    [figH1, lightH] = mbaDisplayConnectome(fg.fibers,figure);
    hold on
    sliceH = feDisplayBrainSlice(anat, [-15 0 0]);
    
    fg = fascicles(3);
    subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)*proportionFibers));
    fg.fibers = fg.fibers(subsampleindx);
    [figH2, lightH2] = mbaDisplayConnectome(fg.fibers,gcf,[.8 .6 .2],'single');
    view(-90,0);
    set(gca,'ylim',[-70 55],'xlim',[-75 0],'zlim',[-25 85])
    delete(lightH)
    delete(lightH2)
    lightH = camlight('right');drawnow;
    saveFig(gcf,fullfile(fig_saveDir,figName))

    % Now plot coronal view:  
    figName = sprintf('left_fascicles_ARC_CST_coronal_%s',lmax{ii});
    view(0,0);
    hold on
    sliceH = feDisplayBrainSlice(anat, [0 -10 0]);
    set(gca,'ylim',[-70 55],'xlim',[-75 0],'zlim',[-25 85])
    delete(lightH)
    lightH = camlight('right');drawnow; 
    saveFig(gcf,fullfile(fig_saveDir,figName))

   % Superior-longitudinal fasciculum and corticospinal tract
    figName = sprintf('right_fascicles_SLF_sagital_%s',lmax{ii});
        
    fg = fascicles(16);
    subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)*proportionFibers));
    fg.fibers = fg.fibers(subsampleindx);

    [figH, lightH] = mbaDisplayConnectome(fg.fibers,figure);
    hold on
    sliceH = feDisplayBrainSlice(anat, [15 0 0]);    
    view(90,0);
    set(gca,'ylim',[-70 55],'xlim',[0 75],'zlim',[-25 85])
    delete(lightH)
    lightH = camlight('right');drawnow;    
    saveFig(gcf,fullfile(fig_saveDir,figName))
  
    % Superior-longitudinal fasciculum and right arcuate
    figName = sprintf('right_fascicles_ARC_SLF_sagital_%s',lmax{ii});
    fg = fascicles(16);
    subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)*proportionFibers));
    fg.fibers = fg.fibers(subsampleindx);

    [figH, lightH] = mbaDisplayConnectome(fg.fibers,figure);
    hold on
    sliceH = feDisplayBrainSlice(anat, [15 0 0]);
      
    fg = fascicles(20);
    subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)*proportionFibers));
    fg.fibers = fg.fibers(subsampleindx);

    [figH, lightH2] = mbaDisplayConnectome(fg.fibers,gcf,[.8 .6 .2],'single');
    view(90,0);  
    set(gca,'ylim',[-70 55],'xlim',[0 75],'zlim',[-25 85])
    delete(lightH)
    delete(lightH2)
    lightH = camlight('right');drawnow; 
    saveFig(gcf,fullfile(fig_saveDir,figName))

    close all

  end
end % Repeated tracking
end % b-value

end % End main function

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

