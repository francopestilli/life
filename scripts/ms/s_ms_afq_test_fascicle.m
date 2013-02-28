function s_ms_afq_test_fascicle
%
% Uses AFQ to segment a connectome and generate 20 major faascicles. 
%
% Then extract one fascicle (arcuate) fasciculum and test its importance
% within the volume take by the fasccile.
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
sdCutoff             = [3 3 3.5 3.5];     % One per conenctome, generally smaller for smaller lmax values
clobber              = [0 0 0 0 1 0]; % Owerwrite all the files.

% DIRECTORY TO LOAD FILES FROM:
% DWI data
dataRootPath  = fullfile('/biac2','wandell6','data','frk','life_dti','FP20120420');
subfolders    = fullfile('150dirs_b1000_1');
baseDir       = fullfile(dataRootPath,subfolders);
dtFile        = fullfile(baseDir,'dt6.mat');
dwiFile       = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
dwiFileRepeat = fullfile(dataRootPath,'raw','0011_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
t1File        = fullfile(dataRootPath,'t1','t1.nii.gz');

% ROIs connectomes
projectDir           = '/azure/scr1/frk/150dirs_b1000_b2000_b4000';
connectSubfolders    = {'life_mrtrix_rep1','life_mrtrix_rep2','life_mrtrix_rep3'};
connectomeFile       = { '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax12_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_prob-500000.pdb'};%, ...

% DIRECTORIES and FILES TO SAVE2:
saveDir              =   fullfile(projectDir,'results');
fe_structuresDir     =  'fe_arcuate_importance';
fas_structuresDir    =  'fas_arcuate_cst_test';
roi_saveDir          =  'roi_arcuate_importance';
arcuateRoiFileName   =  'arcuate_roi_lmax12_prob_importance';  


% FASCICLES:
fascicleNames   = {'rArc', 'lArc','rSlf', 'lSlf'};
fascicleIndices = {20,19,16,15};


% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

for irep = 1:length(connectSubfolders)
  for ii = 1:length(connectomeFile)
    fprintf('[%s] Processing: \n %s \n ======================================== \n\n',mfilename,connectomeFile{ii})
    
    % Name of the whole-brain connectome to load.
    wholeBrainConnectome = fullfile(projectDir,connectSubfolders{1},connectomeFile{1});
    
    % Set up file names
    % Find an identificative name for the connectome that is shortenough:
    cName         = [connectomeFile{ii}(1:57),'_',connectomeFile{ii}(end-16:end-4)];
    
    % Name and path for savign the fe structures
    feSaveDir     = fullfile(saveDir,connectSubfolders{irep},fe_structuresDir);
    fasSaveDir    = fullfile(saveDir,connectSubfolders{irep},fas_structuresDir);
    roiSaveDir    = fullfile(saveDir,connectSubfolders{irep},roi_saveDir);
    
    % Build the full name of the two fascicles FE's structures
    feSaveNameAll     = sprintf('%s_diffModAx%sRd%s_%s',cName,num2str(100*diffusionModelParams(1)), ...
                                                              num2str(100*diffusionModelParams(2)));
    feSaveNameArcuate = sprintf('%s_%s_sd%2.0f',cName,feSaveNameAll,100*sdCutoff(1));
    feSaveNameNOArc   = sprintf('%s_NOT_%s_sd%2.0f',cName,feSaveNameAll,100*sdCutoff(1));
    feSaveNameWITHArc   = sprintf('%s_WITH_%s_sd%2.0f',cName,feSaveNameAll,100*sdCutoff(1));
    
    % Name and path for saving the fascicles
    fasSaveName        = sprintf('fas_%s',feSaveNameAll);
    fasFullPath        = fullfile(fasSaveDir,fasSaveName);
    fasCleanedSaveName = sprintf('%s_cleaned_sd%i',fasSaveName,100*sdCutoff(1));
    fasCleanedFullPath  = fullfile(fasSaveDir,fasCleanedSaveName);
    
    % Build a name for the roi
    arcuateRoiFileName  = sprintf('arcuateROI_%s_sd%2.0f.mat',arcuateRoiFileName,100*sdCutoff(1));
    roiFullPath         = fullfile(roiSaveDir,arcuateRoiFileName); 
 
    % Check if the folders need to be created
    checkFolders({roiFullPath,feSaveNameArcuate,feSaveDir,fasSaveDir,roiSaveDir});

    %% Check if the these fascicles were computed already
    if ~(exist([fasFullPath,'.mat'],'file') == 2) || clobber(1)
      fprintf('[%s] Creating Fascicles: \n%s\n ======================================== \n\n',mfilename,fasFullPath)
      
      % Segment the fiber group using AFQ.
      [fascicles,classification,fg] = feAfqSegment(dtFile, wholeBrainConnectome);
      
      % Save the fascicles and the cleaned whole-brain connectome
      save([fasFullPath,'.mat'],'fascicles','classification','fg','-v7.3');
      clear fg 
      fprintf('[%s] DONE creating fascicles: \n%s\n ======================================== \n\n',mfilename,fasFullPath)
    else
      fprintf('[%s] FOUND Fascicles FILE NOT PROCESSING, Loading it: \n%s\n ======================================== \n\n',mfilename,fasFullPath)
      load([fasFullPath,'.mat'],'fascicles','classification','fg');
    end
    
    %% Check if we have fascicles that were precomputed and cleaned
    if ~(exist([fasCleanedFullPath,'.mat'],'file') == 2) || clobber(2)
      fprintf('[%s] Cleaning fascicles (sd %2.2f): \n%s\n ======================================== \n\n',mfilename,fasFullPath,sdCutoff)
      
      % These are the options for cleaning the outliers from the fascicles.
      % We pass a larger Z-score for the lmax=2, so that the fascicle is
      % larger and has more likelyhood to overlap the CST
      opts.stdCutOff   = sdCutoff(ii);   % Standard deviation fo the 3D gaussian distribution used to
      % represent the fascicle when removing outliers. 3.5 z-scores
      opts.maxLen      = 20;  % Max lenght of fibers to be accepted in cm
      opts.maxNumNodes = 100; % This is used only during the computations does not actually change the nodes
      
      % Remove the outliers from the fascicles
      [fascicles classificationCleaned fg] = feAfqRemoveFascicleOutliers(fascicles, classification,opts);
      
      % Save the fascicles and the cleaned whole-brain connectome
      save([fasCleanedFullPath,'.mat'],'fascicles','classificationCleaned','fg','-v7.3');
      fprintf('[%s] DONE cleaning fascicles: \n%s\n ======================================== \n\n',mfilename,fasFullPath)
    else
      fprintf('[%s] FOUND Cleaned Fascicles NOT PROCESSING, Loading it: \n%s\n ======================================== \n\n',mfilename,fasFullPath)
      load([fasCleanedFullPath,'.mat'],'fascicles','classificationCleaned','fg');
    end
    
    % Extract th two fascicles necessary
    arcuateFG = fascicles(19);% Left arcuate fasciculum
    clear fascicles 

    % Show the fasciles
    %h = figure;feConnectomeDisplay(arcuateFG,h,[.5 .7 .9]); drawnow
      
    %% Build the LiFE model around the volume of the Left Arcuate Fasciculum
    if ~(exist(fullfile(feSaveDir,[feSaveNameArcuate,feSaveNameArcuate,'.mat']),'file') == 2) || clobber(3)
      fprintf('[%s] BUILDING LiFE for arcuate only \n%s\n ======================================== \n\n',mfilename,feSaveNameArcuate)
      fe = feConnectomeInit(dwiFile,dtFile,arcuateFG,feSaveNameArcuate,feSaveDir,dwiFileRepeat, ...
                                    t1File,diffusionModelParams);
      
      % Fit the model with global weights.
      fe    = feSet(fe,'fit',feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsig demeaned'),'sgdnn'));
      
      % Fit the model with voxel-wise weights.
      %fe = feFitModelByVoxel(fe);
      feConnectomeSave(fe,feSaveNameArcuate);
      clear fe
      
      fprintf('[%s] DONE building FE for arcuate only: \n%s\n ======================================== \n\n',mfilename,feSaveNameArcuate)
    else
      fprintf('[%s] FOUND arcuate fe FILE NOT PROCESSING, Loading it: \n%s\n ======================================== \n\n',mfilename,feSaveNameArcuate)
      load(fullfile(feSaveDir,[feSaveNameArcuate,feSaveNameArcuate,'.mat']));
    end
    
    %% Get the brain VOLUME inside which to evaluate the model
    fprintf('[%s] Getting ROI: \n%s\n ======================================== \n\n',mfilename,arcuateRoiFileName)
    if ~exist([roiFullPath],'file') || clobber(4)
      % We extract the ROI of the First fascicle (arcuate) and LMAX value (2, stream for the moment) in ACPC
      % because the fascicles are in ACPC. We save it out as a file. This first ROI will be re-loaded an used
      % for all the rest of the analysis.
      arcuateROI = dtiNewRoi(arcuateRoiFileName, rand(1,3), fefgGet(arcuateFG,'unique image coords')');
      dtiWriteRoi(arcuateROI, roiFullPath, [], 'acpc');
    else
      % Load it the ROI must have been previousaly generated
      arcuateROI = dtiReadRoi(roiFullPath);
    end
    fprintf('[%s] Using ROI: \n%s\n ======================================== \n\n',mfilename,arcuateRoiFileName)
    
    %% Build the LiFE model for the FG WITH the arcuate fasciculum
    if ~(exist(fullfile(feSaveDir,[feSaveNameWITHArc,feSaveNameWITHArc,'.mat']),'file') == 2) || clobber(5)
      % Clip the corticospinal tract to be constrained inside the volume defined
      % by the arcuate.
      tic, fprintf('[%s] Clipping Connectome fibers fibers volume ROI of the arcuate... ',mfilename);
      fg = feClipFibersToVolume(fg,arcuateROI.coords,maxVolDist);
      fprintf('process completed in %2.3fminutes\n',toc/60);
      clear arcuateROI
      
      fprintf('[%s] BUILDING LiFE for connectome without the Arcute \n%s\n ======================================== \n\n',mfilename,feSaveNameWITHArc)
      % Build a new LiFE model with it
      fe   = feConnectomeInit(dwiFile,dtFile,fg,feSaveNameWITHArc,feSaveDir,dwiFileRepeat, ...
        t1File,diffusionModelParams);
      
      % Fit the model with global weights.
      fe    = feSet(fe,'fit',feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsig demeaned'),'sgdnn'));
      
      % Fit the model with voxel-wise weights.
      fe = feFitModelByVoxel(fe);
      feConnectomeSave(fe,feSaveNameWITHArc);
      clear fe
      fprintf('[%s] DONE BUILDING LiFE for the connectome WITH the arcuate fasciculum \n%s\n ======================================== \n\n',mfilename,feSaveNameWITHArc)
    else
      fprintf('FOUND Union (arc/CST) fe FILE NOT PROCESSING, Loading it: \n%s\n ======================================== \n\n',mfilename,feSaveNameWITHArc)
      %load(fullfile(feSaveDir,[feSaveNameArcCST,feSaveNameArcCST,'.mat']));
    end
    
    %% Build the LiFE model for the FG excluding the arcuate fasciculum
    if ~(exist(fullfile(feSaveDir,[feSaveNameNOArc,feSaveNameNOArc,'.mat']),'file') == 2) || clobber(6)
      % Load the whole-brain connectome and the indexing into it.
      load([fasCleanedFullPath,'.mat'],'fg','classificationCleaned');

      % exclude the arcuate fibers:
      % We stored the indices into FG (the whole-brain connectome inside)
      % the 'classification' structure. Now we use tht information to
      % extract the rest of the connectome, without the arcuate fasciculum
      % from the whole-brain conenctome.
      arcuateIndxIntoFG = find(classificationCleaned.index==19);
      notArcuateFG      = fgExtract(fg,arcuateIndxIntoFG,'remove');
      clear classificationClened fg arcuateIndxIntoFG
      
      % Clip the corticospinal tract to be constrained inside the volume defined
      % by the arcuate.
      tic, fprintf('[%s] Clipping Connectome fibers fibers volume ROI of the arcuate... ',mfilename);
      notArcuateFG = feClipFibersToVolume(notArcuateFG,arcuateROI.coords,maxVolDist);
      fprintf('process completed in %2.3fminutes\n',toc/60);
            
      fprintf('[%s] BUILDING LiFE for connectome without the Arcute \n%s\n ======================================== \n\n',mfilename,feSaveNameNOArc)
      % Build a new LiFE model with it
      fe   = feConnectomeInit(dwiFile,dtFile,notArcuateFG,feSaveNameNOArc,feSaveDir,dwiFileRepeat, ...
                               t1File,diffusionModelParams);
      
      % Fit the model with global weights.
      fe    = feSet(fe,'fit',feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsig demeaned'),'sgdnn'));
      
      % Fit the model with voxel-wise weights.
      fe = feFitModelByVoxel(fe);
      feConnectomeSave(fe,feSaveNameNOArc);
      clear fe
      fprintf('[%s] DONE BUILDING LiFE for the connectome excluding th arcuate fasciculum \n%s\n ======================================== \n\n',mfilename,feSaveNameNOArc)
    else
      fprintf('FOUND Union (arc/CST) fe FILE NOT PROCESSING, Loading it: \n%s\n ======================================== \n\n',mfilename,feSaveNameNOArc)
      %load(fullfile(feSaveDir,[feSaveNameArcCST,feSaveNameArcCST,'.mat']));
    end
    
  end
end

end % End function


%-----------------------%
function fold = checkFolders(foldersToCheck)
% Make sure the folders exist otherwise create them:
for ff = 1:length(foldersToCheck)
  [fold,~,~] = fileparts(foldersToCheck{ff});
  if ~( exist(fold,'dir') == 7)
    mkdir(fold)
  end
end
end




