function s_ms_afq_test_fascicle_prediction_volumes
%
% Generates a volume of predicted signal for one of the fascicles geenrated with afq. 
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
maxVolDist           = 2;           % Max distance in mm from the ROI edges.
sdCutoff             = [2.5];     % One per conenctome, generally smaller for smaller lmax values

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
connectSubfolders    = {'life_mrtrix_rep1'};%,'life_mrtrix_rep2','life_mrtrix_rep3'};
connectomeFile       = { '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax12_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_prob-500000.pdb'};

% DIRECTORIES and FILES TO SAVE:
saveDir              =   fullfile(projectDir,'results');
fe_structuresDir     =  'fe_arcuate_importance';
fas_structuresDir    =  'fas_arcuate_cst_test';
roi_saveDir          =  'roi_arcuate_importance';
arcuateRoiFileName   =  'arcuate_roi_lmax12_prob_importance';  


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
    feSaveNameArcuate = sprintf('arcuateFE_%s_sd%2.0f',feSaveNameAll,100*sdCutoff(1));
    feSaveNameNOArc   = sprintf('arcuateFE_NOT_%s_sd%2.0f',feSaveNameAll,100*sdCutoff(1));
    feSaveNameWITHArc = sprintf('arcuateFE_WITH_%s_sd%2.0f',feSaveNameAll,100*sdCutoff(1));
    
    % Name and path for saving the fascicles
    fasSaveName        = sprintf('fas_%s',feSaveNameAll);
    fasFullPath        = fullfile(fasSaveDir,fasSaveName);
    fasCleanedSaveName = sprintf('%s_cleaned_sd%i',fasSaveName,100*sdCutoff(1));
    fasCleanedFullPath = fullfile(fasSaveDir,fasCleanedSaveName);
    
    % Build a name for the roi
    arcuateRoiFileName  = sprintf('arcuateROI_%s_sd%2.0f.mat',arcuateRoiFileName,100*sdCutoff(1));
    roiFullPath         = fullfile(roiSaveDir,arcuateRoiFileName); 
       
    %% Build the LiFE model around the volume of the Left Arcuate Fasciculum
    if ~(exist(fullfile(feSaveDir,[feSaveNameArcuate,feSaveNameArcuate,'.mat']),'file') == 2) 
      fprintf('[%s] CANNOT FIND the LiFE model for the arcuate only \n%s\n ======================================== \n\n',mfilename,feSaveNameArcuate)
      keyboard
    else
      fprintf('[%s] FOUND fe FILE ofr acruate ONLY (no connectome), Loading it: \n%s\n ======================================== \n\n',mfilename,feSaveNameArcuate)
      feONLY = load(fullfile(feSaveDir,[feSaveNameArcuate,feSaveNameArcuate,'.mat']),'fe');
      feONLY = feONLY.fe;
    end
    

    %% Build the LiFE model for the FG WITH the arcuate fasciculum
    if ~(exist(fullfile(feSaveDir,[feSaveNameWITHArc,feSaveNameWITHArc,'.mat']),'file') == 2)
      fprintf('[%s] CANNOT FIND the LiFE model for the connectome WITH the arcuate fasciculum \n%s\n ======================================== \n\n',mfilename,feSaveNameWITHArc)
      keyboard
    else
      fprintf('[%s] FOUND fe FILE with Arcuate plus conenctome Loading it: \n%s\n ======================================== \n\n',mfilename,feSaveNameWITHArc)
      feWITH = load(fullfile(feSaveDir,[feSaveNameWITHArc,feSaveNameWITHArc,'.mat']),'fe');
      feWITH = feWITH.fe;
    end
    
    %% Build the LiFE model for the FG excluding the arcuate fasciculum
    if ~(exist(fullfile(feSaveDir,[feSaveNameNOArc,feSaveNameNOArc,'.mat']),'file') == 2)
      fprintf('[%s] CANNOT FIND the LiFE model for the connectome excluding th arcuate fasciculum \n%s\n ======================================== \n\n',mfilename,feSaveNameNOArc)
      keyboard
    else
      fprintf('[%s] FOUND fe FILE without arcuate, Loading it: \n%s\n ======================================== \n\n',mfilename,feSaveNameNOArc)
      feNOT = load(fullfile(feSaveDir,[feSaveNameNOArc,feSaveNameNOArc,'.mat']),'fe');
      feNOT = feNOT.fe;
    end    %% Get the brain VOLUME inside which to evaluate the model
    
    if ~exist([roiFullPath],'file')
      fprintf('[%s] CANNOT FIND the ROI file \n%s\n ======================================== \n\n',mfilename,arcuateRoiFileName)
      keyboard
    else
      fprintf('[%s] Loading ROI: \n%s\n ======================================== \n\n',mfilename,arcuateRoiFileName)
      % Load it the ROI must have been previousaly generated
      arcuateROI = dtiReadRoi(roiFullPath);
    end
    
    % make some plots
    statType = {'vox r2','vox r2 voxel wise','vox rmse','vox rmse voxel wise','vox rmse ratio','vox rmse ratio voxel wise'};
    for ss = 1:length(statType)
      val = [median(feGetRep(feONLY,statType{ss})), median(feGetRep(feNOT,statType{ss})), median(feGetRep(feWITH,statType{ss}))];
      
      mrvNewGraphWin(upper(statType{ss}));
      bar(val,'k'); 
      ylabel(sprintf('%s\n(Cross validated)',statType{ss}),'fontsize',20)
      set(gca,'xticklabel',{'Only Arcuate','Without Arcuate','With Arcuate'},...
        'box','off','tickdir','out','ticklen',[.01 .05],'FontSize',20);
    end
    
   keyboard 
  end
end

end % End function
