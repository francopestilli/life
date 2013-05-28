function s_ms_afq_test_arcuate_cst_NEW
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
% Written by Franco Pestilli (c) Stanford Vista Team 2013

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
clobber              = [0 0 0 0 0]; % Owerwrite all the files.
plotFG  =0; % Shows he fiber group for the arcuate.
% DIRECTORY TO LOAD FILES FROM:
% DWI data
projectDir  = '/azure/scr1/frk/150dirs_b1000_b2000_b4000';
connectSubfolders    = {'life_mrtrix_rep1','life_mrtrix_rep2','life_mrtrix_rep3'};
% Make some plots
colors = {[.85 .5 .65],[.6 .85 .4],[.45 .6 .95]};
if notDefined('fig_saveDir'), fig_saveDir = fullfile('/home/frk/Dropbox','arcuate_cst_test_figures');end

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
  connectomeFile= {'0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax2_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb'};
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
    if plotFG
        fg = fascicles(19);
        subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)/12));
        fg.fibers = fg.fibers(subsampleindx);
        feConnectomeDisplay(fg,figure, [.3 .7 .9]);
        view(-80,0);camlight('right');drawnow;
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
    
    % Extract th two fascicles necessary
    arcuateFG = fascicles(19);% Left arcuate fasciculum
    cstFG     = fascicles(3); % Left corticospinal tract
    clear fascicles
    
    if plotFG
        fg = arcuateFG;
        subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)/6));
        fg.fibers = fg.fibers(subsampleindx);
        feConnectomeDisplay(fg,figure, [.3 .9 .4]);
        view(-80,0);camlight('right');drawnow;
    end
    
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

      feArc.fe = feConnectomeInit(dwiFile,dtFile,arcuateFG,feSaveNameArcuate,feSaveDir,dwiFileRepeat, ...
                                    t1File,diffusionModelParams);
           
      fg = feGet(feArc.fe,'fg img'); 
      subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)/6));
      fg.fibers = fg.fibers(subsampleindx);
      feConnectomeDisplay(fg,figure, [.9 .7 .5]);
      view(-80,0);camlight('right');drawnow;
      
      % Fit the model with global weights.
      feArc.fe = feSet(feArc.fe,'fit',feFitModel(feGet(feArc.fe,'Mfiber'),feGet(feArc.fe,'dsig demeaned'),'sgdnn'));
      
      % Fit the model with voxel-wise weights.
%      fe = feFitModelByVoxel(fe);
       feConnectomeSave(feArc.fe,feSaveNameArcuate);
      clear fg
      
      fprintf('[%s] DONE  for arcuate ONLY \n',mfilename)
    else
      fprintf('[%s] FOUND arcuate fe FILE NOT PROCESSING, Loading it\n',mfilename)
      feArc = load(fullfile(feSaveDir,[feSaveNameArcuate,feSaveNameArcuate,'.mat']));
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
      clear arcuateFG 
      
      fprintf('[%s] BUILDING LiFE for Union (arc/CST)...',mfilename)
      % Build a new LiFE model with it
      feARC_CST.fe   = feConnectomeInit(dwiFile,dtFile,unionFG,feSaveNameArcCST,feSaveDir,dwiFileRepeat, ...
                               t1File,diffusionModelParams);
      clear unionFG
               
      fg = feGet(feARC_CST.fe,'fg img');
      subsampleindx = randsample(length(fg.fibers),floor(length(fg.fibers)/6));
      fg.fibers = fg.fibers(subsampleindx);
      feConnectomeDisplay(fg,figure, [.9 .8 .4]);
      view(-80,0);camlight('right');drawnow;
      
      % Fit the model with global weights.
      feARC_CST.fe    = feSet(feARC_CST.fe,'fit',feFitModel(feGet(feARC_CST.fe,'Mfiber'),feGet(feARC_CST.fe,'dsig demeaned'),'sgdnn'));
      
      % Fit the model with voxel-wise weights.
     % fe = feFitModelByVoxel(fe);
      feConnectomeSave(feARC_CST.fe,feSaveNameArcCST);
      clear fg
      fprintf('[%s] DONE BUILDING LiFE for Union (arc/CST) \n',mfilename)
    else
      fprintf('[%s] FOUND Union (arc/CST) fe FILE NOT PROCESSING, Loading it\n',mfilename)
      feARC_CST = load(fullfile(feSaveDir,[feSaveNameArcCST,feSaveNameArcCST,'.mat']));
    end
    
    %% Test the improvement in Rrmse only in the voxels shared by Arcuate and CST 
    % Reduce the voxels of the two FE structures to the union of the voxels of
    % the arcuate and the CST
       
    % (1) Make sure the FE with and with CST have the same voxels:
    commonVoxels = ismember(feGet(feARC_CST.fe,'roi coords'),feGet(feArc.fe,'roi coords'),  'rows');
    [feARC_CST.fe, indicesFibersKept]    = feConnectomeReduceVoxels(feARC_CST.fe, find(commonVoxels));

    % (2) find the voxels of the corticospinal tract in the FE structure of
    % the arcutae alone. Reduce the acruate only to those voxels, the
    % intersection between arcuate and cst.
    commonVoxels = ismember(feGet(feArc.fe,'roi coords'), ...
        fefgGet(dtiXformFiberCoords(cstFG, feARC_CST.fe.life.xform.acpc2img,'img'), ...
        'unique image coords'), 'rows');
    feWithoutFas = feConnectomeReduceVoxels(feArc.fe,find(commonVoxels));
   
    % Store the rmse ad Rrmse
    WITHOUT.rrmse(irep) = median(feGetRep(feWithoutFas,'vox rmse ratio'));
    WITHOUT.rmse(irep) = median(feGetRep(feWithoutFas,'vox rmse'));
    
    WITHOUT.rrmseall{irep} = (feGetRep(feWithoutFas,'vox rmse ratio'));
    WITHOUT.rmseall{irep} = (feGetRep(feWithoutFas,'vox rmse'));
    
    commonVoxels = ismember(feGet(feARC_CST.fe,'roi coords'), ...
        fefgGet(dtiXformFiberCoords(cstFG, ...
        feARC_CST.fe.life.xform.acpc2img,'img'),'unique image coords'),  'rows');
    [feWithFas, indicesFibersKept]    = feConnectomeReduceVoxels(feARC_CST.fe, find(commonVoxels));
    
    % Store the rmse ad Rrmse
    WITH.rrmse(irep) = median(feGetRep(feWithFas,'vox rmse ratio'));
    WITH.rmse(irep) = median(feGetRep(feWithFas,'vox rmse'));
     
    WITH.rrmseall{irep} = (feGetRep(feWithFas,'vox rmse ratio'));
    WITH.rmseall{irep} = (feGetRep(feWithFas,'vox rmse'));

  end
end % Repeated tracking
keyboard

% Compute a test of the diference in rmse
% (1) Get the differece in rmse observed empiriclly
EmpiricalDiff = WITHOUT.rmse(1) - WITH.rmse(1);

% (2) Compute the Null distribution by:
% (2.1) Combine all the rmse from both WITH and WITHOUT.
% (2.2) Compute 10,000 distributions of rmse one for WITH one for WITHOUT
%       with the same numerosity of the initial samples
% (2.3) Compute the difference between the medians of these 10,000
%       distributions.
NullSet     = [WITHOUT.rmseall{1} WITH.rmseall{1}];
sizeWith    = length(WITH.rmseall{1});
sizeWithout = length(WITHOUT.rmseall{1});

nboots = 100000;
nullDistribution = nan(nboots,1);
parfor ibt = 1:nboots
    % Y = RANDSAMPLE(POPULATION,K)
    BootWith    = median(randsample(NullSet,sizeWith));
    BootWithout = median(randsample(NullSet,sizeWithout));   
    nullDistribution(ibt) = BootWithout - BootWith;
end

% Plot the null distribution and the empirical difference
figName = sprintf('Test_Arcuate_CST_BOOT_test_rmse_%s',cName);
fh = mrvNewGraphWin(figName);
[y,x] = hist(nullDistribution,100);
y = y./sum(y);
bar(x,y,'k')
hold on
plot([EmpiricalDiff,EmpiricalDiff],[0 max(y)],'r-','linewidth',2)
set(gca,'tickdir','out','box','off', 'ylim',[0 max(y)],'FontSize',16)
ylabel('Likelihood')
xlabel('Difference in rmse')

% (3) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
if max(nullDistribution)<EmpiricalDiff
    p = 100*1/nboots;
else
    p = sum(nullDistribution(sort(nullDistribution)>EmpiricalDiff));
end
title(sprintf('The probability of obtaining the difference by chance is less then %2.3f%%',p), ...
    'FontSize',16)
saveFig(fh,fullfile(fig_saveDir,figName))


% MAke aplot across repeated tracking
% Now make averages and std or the rmse, Rrmse, r2 values for plotting
WITH.rmsem   = mean(WITH.rmse);
WITH.rmsesd  = [WITH.rmsem-std(WITH.rmse); WITH.rmsem+std(WITH.rmse)];
WITH.rrmsem  = mean(WITH.rrmse);
WITH.rrmsesd = [WITH.rrmsem-std(WITH.rrmse); WITH.rrmsem+std(WITH.rrmse)];

WITHOUT.rmsem   = mean(WITHOUT.rmse);
WITHOUT.rmsesd  = [WITHOUT.rmsem-std(WITHOUT.rmse); WITHOUT.rmsem+std(WITHOUT.rmse)];
WITHOUT.rrmsem  = mean(WITHOUT.rrmse);
WITHOUT.rrmsesd = [WITHOUT.rrmsem-std(WITHOUT.rrmse); WITHOUT.rrmsem+std(WITHOUT.rrmse)];

% Make a plot of the R-squared
figName = sprintf('Test_Arcuate_CST_Addition_rmse_%s',cName);
fh = mrvNewGraphWin(figName);
bar([WITH.rmsem,WITHOUT.rmsem],'FaceColor',colors{1})
hold on
plot([1 1; 2 2]',[WITH.rmsesd,WITHOUT.rmsesd],'r-','linewidth',16)
ylabel('rmse')
set(gca,'xticklabel',{'Arcuate and CST','Arcuate alone'},'tickdir','out','box','off', ...
    'ylim',[20 50],'FontSize',16)
saveFig(fh,fullfile(fig_saveDir,figName))

% Make a plot of the R-squared
figName = sprintf('Test_Arcuate_CST_Addition_%s',cName);
fh = mrvNewGraphWin(figName);
bar([WITH.rrmsem,WITHOUT.rrmsem],'FaceColor',colors{1})
hold on
plot([1 1; 2 2]',[WITH.rrmsesd,WITHOUT.rrmsesd],'r-','linewidth',16)
plot([0 3],[1 1],'k-')
ylabel('R_{rmse}')
set(gca,'xticklabel',{'Arcuate and CST','Arcuate alone'},'tickdir','out','box','off', ...
    'ylim',[0.75 1.5],'FontSize',16)
saveFig(fh,fullfile(fig_saveDir,figName))

% Make a scatter plto:
figName = sprintf('Test_Arcuate_CST_Addition_rmse_SCATTER_%s',cName);
fh = mrvNewGraphWin(figName);
hold on
plot([0 80],[0 80],'k-',[median(WITHOUT.rmse),median(WITHOUT.rmse)],[0 80],'k--',[0 80],[median(WITH.rmse),median(WITH.rmse)],'k--')

for ii = 1:length(WITHOUT.rmseall)
plot(WITHOUT.rmseall{ii},WITH.rmseall{ii},'o','color',colors{ii},'markerfacecolor',colors{ii});
end
axis square
set(gca,'tickdir','out','box','off','xlim',[0 80],'ylim',[0 80],'FontSize',16)
ylabel('rmse WITH cortico-spinal tract')
xlabel('rmse WITHOUT cortico-spinal tract')
saveFig(fh,fullfile(fig_saveDir,figName))

% Make a scatter plto:
figName = sprintf('Test_Arcuate_CST_Addition_SCATTER_%s',cName);
fh = mrvNewGraphWin(figName);
plot([1 1],[.5 3],'k--',[.5 3],[1 1],'k--',[.5 3],[.5 3],'k-')
hold on
for  ii = 1:length(WITHOUT.rmseall)
plot(WITHOUT.rrmseall{ii},WITH.rrmseall{ii},'ro','color',colors{ii},'markerfacecolor',colors{ii})
end
axis square
set(gca,'tickdir','out','box','off','ylim',[0.5 3],'FontSize',16)
ylabel('R_{rmse} WITH cortico-spinal tract')
xlabel('R_{rmse} WITHOUT cortico-spinal tract')
saveFig(fh,fullfile(fig_saveDir,figName))

keyboard
end

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

eval( sprintf('print(%s,  ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'');', num2str(h),figName));

end

