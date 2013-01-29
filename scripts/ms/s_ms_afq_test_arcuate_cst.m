function s_ms_afq_test_arcuate_cst
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
for bb = 1:length(bval)
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
      
      % Segment the fiber group using AFQ.
      [fascicles,classification,fg] = feAfqSegment(dtFile, wholeBrainConnectome);
      
      % Save the fascicles and the cleaned whole-brain connectome
      save([fasFullPath,'.mat'],'fascicles','classification','fg','-v7.3');
      fprintf('[%s] DONE creating fascicles \n',mfilename)
    else
      fprintf('[%s] FOUND Fascicles FILE NOT PROCESSING, Loading it\n',mfilename)
      load([fasFullPath,'.mat'],'fascicles','classification');
    end
    
    fg = fascicles(19);
    subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)/12));
    fg.fibers = fg.fibers(subsampleindx);
    feConnectomeDisplay(fg,figure, [.3 .7 .9]);
    view(-80,0);camlight('right');drawnow;
    
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
      [fascicles classificationCleaned] = feAfqRemoveFascicleOutliers(fascicles, classification,opts);
      
      % Save the fascicles and the cleaned whole-brain connectome
      save([fasCleanedFullPath,'.mat'],'fascicles','classificationCleaned','-v7.3');
      fprintf('[%s] DONE cleaning fascicles \n',mfilename)
    else
      fprintf('[%s] FOUND Cleaned Fascicles NOT PROCESSING, Loading it\n',mfilename)
      load([fasCleanedFullPath,'.mat'],'fascicles','classificationCleaned');
    end
    
    % Extract th two fascicles necessary
    arcuateFG = fascicles(19);% Left arcuate fasciculum
    cstFG     = fascicles(3); % Left corticospinal tract
    clear fascicles 
      
    fg = arcuateFG;
    subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)/6));
    fg.fibers = fg.fibers(subsampleindx);
    feConnectomeDisplay(fg,figure, [.3 .9 .4]);
    view(-80,0);camlight('right');drawnow;

    %% Get the brain VOLUME inside qhich to evaluate the model
    fprintf('[%s] Loading ROI: \n',mfilename)
    arcuateROI = dtiReadRoi(roiFullPath);

    %% Build the LiFE model around the volume of the Left Arcuate Fasciculum
    if ~(exist(fullfile(feSaveDir,[feSaveNameArcuate,feSaveNameArcuate,'.mat']),'file') == 2) || clobber(3)
      fprintf('[%s] BUILDING LiFE for arcuate only \n',mfilename)
      % Clip the corticospinal tract to be constrained inside the volume defined
      % by the arcuate.
      tic, fprintf('[%s] Clipping Arcuate fibers volume ROI of the arcuate... ',mfilename);
      arcuateFG = feClipFibersToVolume(arcuateFG,arcuateROI.coords,maxVolDist);
      fprintf('process completed in %2.3fminutes\n',toc/60);

      fe = feConnectomeInit(dwiFile,dtFile,arcuateFG,feSaveNameArcuate,feSaveDir,dwiFileRepeat, ...
                                    t1File,diffusionModelParams);
                                  
      fg = feGet(fe,'fg img'); 
      subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)/6));
      fg.fibers = fg.fibers(subsampleindx);
      feConnectomeDisplay(fg,figure, [.9 .7 .5]);
      view(-80,0);camlight('right');drawnow;
      
      % Fit the model with global weights.
      fe = feSet(fe,'fit',feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsig demeaned'),'sgdnn'));
      
      % Fit the model with voxel-wise weights.
      fe = feFitModelByVoxel(fe);
      feConnectomeSave(fe,feSaveNameArcuate);
      clear fe fg
      
      fprintf('[%s] DONE  for arcuate ONLY \n',mfilename)
    else
      fprintf('[%s] FOUND arcuate fe FILE NOT PROCESSING, Loading it\n',mfilename)
      %load(fullfile(feSaveDir,[feSaveNameArcuate,feSaveNameArcuate,'.mat']));
    end
    

    %% Build the LiFE model for the union fascicle (arc+CST)
    if ~(exist(fullfile(feSaveDir,[feSaveNameArcCST,feSaveNameArcCST,'.mat']),'file') == 2) || clobber(5)
      % Clip the corticospinal tract to be constrained inside the volume defined
      % by the arcuate.
      tic, fprintf('[%s] Clipping CST fibers volume ROI of the arcuate... ',mfilename);
      cstFG = feClipFibersToVolume(cstFG,arcuateROI.coords,maxVolDist);
      fprintf('process completed in %2.3fminutes\n',toc/60);
      
      % Make a Union between the Left Arcuate and Cortico spinal tract.
      unionFG = fgMerge(arcuateFG,cstFG,'Union of Left Arcuate and CST');
      clear arcuateFG cstFG
      
      fprintf('[%s] BUILDING LiFE for Union (arc/CST)...',mfilename)
      % Build a new LiFE model with it
      fe   = feConnectomeInit(dwiFile,dtFile,unionFG,feSaveNameArcCST,feSaveDir,dwiFileRepeat, ...
                               t1File,diffusionModelParams);
      clear unionFG
               
      fg = feGet(fe,'fg img');
      subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)/6));
      fg.fibers = fg.fibers(subsampleindx);
      feConnectomeDisplay(fg,figure, [.9 .8 .4]);
      view(-80,0);camlight('right');drawnow;
      
      % Fit the model with global weights.
      fe    = feSet(fe,'fit',feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsig demeaned'),'sgdnn'));
      
      % Fit the model with voxel-wise weights.
      fe = feFitModelByVoxel(fe);
      feConnectomeSave(fe,feSaveNameArcCST);
      clear fe fg
      fprintf('[%s] DONE BUILDING LiFE for Union (arc/CST) \n',mfilename)
    else
      fprintf('[%s] FOUND Union (arc/CST) fe FILE NOT PROCESSING, Loading it\n',mfilename)
      %load(fullfile(feSaveDir,[feSaveNameArcCST,feSaveNameArcCST,'.mat']));
    end
    
  end
end
end
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




