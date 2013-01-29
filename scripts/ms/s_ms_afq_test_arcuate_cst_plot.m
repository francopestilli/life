function s_ms_afq_test_arcuate_cst_plot
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
sdCutoff             = [3.32 4.7];     % One per conenctome, generally smaller for smaller lmax values

% Recompute the results (1) or load them from disk (0)
recomputeResults = 1;

% MOVIE INFO
makeMovies   = 1;
numRotations = 38; % 180/38*numAnglesPerRotations (which iss et to 5 in feMakeMovie)

% DIRECTORY TO LOAD FILES FROM:
% DWI data
projectDir        = '/azure/scr1/frk/150dirs_b1000_b2000_b4000';
connectSubfolders = {'life_mrtrix_rep1','life_mrtrix_rep2','life_mrtrix_rep3'};

% Type of statistcs to compute, we compute 2 each for the arcuate alone and
% the arcuate plus the CST
stats = {'vox rmse ratio','vox rmse','vox rmse ratio voxel wise','vox rmse voxel wise'};

% DIRECTORIES and FILES TO SAVE2:
saveDir              =   fullfile(projectDir,'results');
fe_structuresDir     =  'fe_arcuate_cst_test';
fas_structuresDir    =  'fas_arcuate_cst_test';
roi_saveDir          =  'roi_arcuate';
arcuateRoiFileName   =  'arcuate_roi_lmax12_prob';
plots_saveDir        = fullfile('/home/frk/Dropbox','arcuate_cst_plots');
movies_saveDir        = fullfile('/home/frk/Dropbox','arcuate_cst_movies');

if recomputeResults
  % Handling parallel processing
  poolwasopen=1; % if a matlabpool was open already we do not open nor close one
  if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end
  
  % Select the dwi file
  bval  = [1000 2000 4000];
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
    
    for i_lmax = 1:length(connectomeFile)
      fprintf('\n\n[%s] NEW CONNECTOME ###############################################\n',mfilename);
      for irep = 1:length(connectSubfolders)
        
        % Name of the whole-brain connectome to load.
        wholeBrainConnectome = fullfile(projectDir,connectSubfolders{irep},connectomeFile{i_lmax});
        
        % Set up file names
        % Find an identificative name for the connectome that is shortenough:
        cName = [connectomeFile{i_lmax}(1:57),'_',connectomeFile{i_lmax}(end-16:end-4)];
        
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
        fasCleanedSaveName = sprintf('%s_cleaned_sd%i',fasSaveName,100*sdCutoff(i_lmax));
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
        
        % Make movies of the Fascicles
        if makeMovies
          % Check if we have fascicles that were precomputed and cleaned
          if ~(exist([fasCleanedFullPath,'.mat'],'file') == 2)
            fprintf('[%s] Cannot find cleaned fascicles \n',mfilename)
            keyboard
          else
            fprintf('[%s] FOUND Cleaned Fascicles NOT PROCESSING, Loading it\n',mfilename)
            load([fasCleanedFullPath,'.mat'],'fascicles');
          end
          makeFasciclesMovies(fascicles(3),fascicles(19),numRotations,movies_saveDir,fasCleanedSaveName)
          clear fascicles
        end
        
        %% Get the brain VOLUME inside qhich to evaluate the model
        %     fprintf('[%s] Loading ROI: \n',mfilename)
        %     arcuateROI = dtiReadRoi(roiFullPath);
        %     arcuateROI.coords = arcuateROI.coords';
        %     AFQ_RenderRoi(arcuateROI, [.9 .8 .4], 'isosurf', 'wire');
        %     set(gcf,'color','k');set(gca,'color','k')
        %     view(-80,0);  camlight right
        
        %% Build the LiFE model around the volume of the Left Arcuate Fasciculum
        if ~(exist(fullfile(feSaveDir,[feSaveNameArcuate,feSaveNameArcuate,'.mat']),'file') == 2)
          fprintf('[%s] Cannot find FE for the arcuate fasciculum ONLY \n',mfilename)
          keyboard
        else
          fprintf('[%s] FOUND arcuate fe FILE Loading it\n',mfilename)
          load(fullfile(feSaveDir,[feSaveNameArcuate,feSaveNameArcuate,'.mat']));
          feArc = fe; clear fe
        end
        
        
        %% Build the LiFE model for the union fascicle (arc+CST)
        if ~(exist(fullfile(feSaveDir,[feSaveNameArcCST,feSaveNameArcCST,'.mat']),'file') == 2)
          fprintf('[%s] Cannot find the FE structure for the Union (arc/CST) \n',mfilename)
          keyboard
        else
          fprintf('[%s] FOUND Union (arc/CST) FE file, Loading it\n',mfilename)
          load(fullfile(feSaveDir,[feSaveNameArcCST,feSaveNameArcCST,'.mat']));
          feUnion = fe; clear fe
        end
        
        %% Compute some statistcs
        for is = 1:length(stats)
          % Arcuate fascicle only
          a(bb,i_lmax,is,irep)   = median(feGetRep(feArc,  stats{is}));
          
          % Union Fascicle
          u(bb,i_lmax,is,irep)   = median(feGetRep(feUnion,stats{is}));
        end
        
      end
    end   
  end
  
  % Save the results colelcted so far
  save(fullfile(plots_saveDir,sprintf('%s_results.mat',mfilename)),'a','u')
else
  load(fullfile(plots_saveDir,sprintf('%s_results.mat',mfilename)),'a','u')
  
end

% Make plots
for ii = 1:length(stats)
  makeBarPlotAcrossLmax(a,u, plots_saveDir,stats,ii)
  makeBarPlotAcrossBvals(a,u,plots_saveDir,stats,ii)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeBarPlotAcrossLmax(a,u,plots_saveDir,stats,ii)
% Make a plot across bval
colors = {[.25 .4 .45],[.35 .4 .55],[.5 .4 .65]};
xa     = {[1 1 2 2],[4 4 5 5]};
xu     = {[2 2 3 3],[5 5 6 6]};
bval   =  [1000,2000,4000];

switch stats{ii}
  case {'vox rmse ratio'}
    ymin = 0.9;ymax = 1.5;
    ticks = [ymin 1 ymax];
    
  case {'vox rmse ratio voxel wise'}
    ymin = 0.85;ymax = 1.15;
    ticks = [ymin 1 ymax];
    
  case {'vox rmse'}
    ymin = 25;ymax = 55;
    ticks = [ymin (ymax+ymin)/2 ymax];
    
  case {'vox rmse voxel wise'}
    ymin = 25;ymax = 45;
    ticks = [ymin (ymax+ymin)/2 ymax];
    
  otherwise
    keyboard
end

for bb = 1:3
  figName = sprintf('%s_%s_bval%i_lmax_12_2',mfilename,stats{ii}(~isspace(stats{ii})),bval(bb));
  h = mrvNewGraphWin(figName);clf; hold on
  plot([0 9],[1 1],'k--')
  set(gca,'box','off','tickDir','out','ticklen',[0.014 0], ...
    'ylim',[ymin ymax],...
    'ytick',ticks, ...
    'xlim',[0 7],...
    'xtick',[2 5], ...
    'xticklabel',{'Lmax=12','Lmax=2'})
  
  for i_lmax = 1:2
    patch([xa{i_lmax}], [ymin mean(squeeze(a(bb,i_lmax,ii,:))) mean(squeeze(a(bb,i_lmax,ii,:))) ymin], ...
      'k','FaceColor',colors{i_lmax},...
      'EdgeColor','w');
    plot([mean(xa{i_lmax}) mean(xa{i_lmax})],[mean(squeeze(a(bb,i_lmax,ii,:)))-2*std(squeeze(a(bb,i_lmax,ii,:))), ...
          mean(squeeze(a(bb,i_lmax,ii,:)))+2*std(squeeze(a(bb,i_lmax,ii,:)))], ...
      'r-','linewidth',8);
    
    patch([xu{i_lmax}], [ymin mean(squeeze(u(bb,i_lmax,ii,:))) mean(squeeze(u(bb,i_lmax,ii,:))) ymin], ...
      'k','FaceColor',colors{i_lmax},...
      'EdgeColor','w');
    plot([mean(xu{i_lmax}) mean(xu{i_lmax})],[mean(squeeze(u(bb,i_lmax,ii,:)))-2*std(squeeze(u(bb,i_lmax,ii,:))), ...
      mean(squeeze(u(bb,i_lmax,ii,:)))+2*std(squeeze(u(bb,i_lmax,ii,:)))], ...
      'r-','linewidth',8);
    
  end
  % Save the figure
  saveFig(h,fullfile(plots_saveDir,figName))
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeBarPlotAcrossBvals(a,u,plots_saveDir,stats,ii)
% Make a plot across bval
colors = {[.25 .4 .45],[.35 .4 .55],[.5 .4 .65]};
xa     = {[1 1 2 2],[2 2 3 3],[3 3 4 4]};
xu     = {[5 5 6 6],[6 6 7 7],[7 7 8 8]};
lmax = [12,2];

switch stats{ii}
  case {'vox rmse ratio'}
    ymin = 0.9;ymax = 1.5;
    ticks = [ymin 1 ymax];
    
  case {'vox rmse ratio voxel wise'}
    ymin = 0.85;ymax = 1.15;
    ticks = [ymin 1 ymax];
    
  case {'vox rmse'}
    ymin = 25;ymax = 55;
    ticks = [ymin (ymax+ymin)/2 ymax];
    
  case {'vox rmse voxel wise'}
    ymin = 25;ymax = 45;
    ticks = [ymin (ymax+ymin)/2 ymax];
    
  otherwise
    keyboard
end

for i_lmax = 1:2
  figName = sprintf('%s_%s_lmax%i_bval_1000_2000_4000',mfilename,stats{ii}(~isspace(stats{ii})),lmax(i_lmax));
  h = mrvNewGraphWin(figName);clf
  hold on
  plot([0 9],[1 1],'k--')
  set(gca,'box','off','tickDir','out','ticklen',[0.014 0], ...
    'ylim',[ymin ymax],...
    'ytick',ticks, ...
    'xlim',[0 9],...
    'xtick',[2.5 6.5], ...
    'xticklabel',{'Arcuate','Arcuate + CST'})
  
  for bb = 1:3
    patch([xa{bb}], [ymin mean(squeeze(a(bb,i_lmax,ii,:))) mean(squeeze(a(bb,i_lmax,ii,:))) ymin], ...
      'k','FaceColor',colors{bb},...
      'EdgeColor','w');
    plot([mean(xa{bb}) mean(xa{bb})],[mean(squeeze(a(bb,i_lmax,ii,:)))-2*std(squeeze(a(bb,i_lmax,ii,:))), ...
          mean(squeeze(a(bb,i_lmax,ii,:)))+2*std(squeeze(a(bb,i_lmax,ii,:)))], ...
      'r-','linewidth',8);
    
    patch([xu{bb}], [ymin mean(squeeze(u(bb,i_lmax,ii,:))) mean(squeeze(u(bb,i_lmax,ii,:))) ymin], ...
      'k','FaceColor',colors{bb},...
      'EdgeColor','w');
    plot([mean(xu{bb}) mean(xu{bb})],[mean(squeeze(u(bb,i_lmax,ii,:)))-2*std(squeeze(u(bb,i_lmax,ii,:))), ...
      mean(squeeze(u(bb,i_lmax,ii,:)))+2*std(squeeze(u(bb,i_lmax,ii,:)))], ...
      'r-','linewidth',8);
    
  end
  % Save the figure
  saveFig(h,fullfile(plots_saveDir,figName))
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)

printCommand = sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName);
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

% do the printing here:
eval(printCommand);

% save a wiki-compatible file
eval(sprintf('print(%s, ''-painters'',''-dpng'', ''-noui'', ''%s'')', num2str(h),figName));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeFasciclesMovies(cstFG,arcuateFG,numRotations,fasSaveDir,fasCleanedSaveName)

% Show the fasciles
h = figure;
feConnectomeDisplay(arcuateFG,h,[.5 .7 .9]); hold on
view(-80,0); camlight right; drawnow
feMakeMovie(h,numRotations,fullfile(fasSaveDir,sprintf('%s_arcuate_%s',mfilename,fasCleanedSaveName)))
close(h)

% Show the fasciles
h = figure;
feConnectomeDisplay(cstFG,h,[.9 .7 .5]); hold on
view(-80,0); camlight right; drawnow
feMakeMovie(h,numRotations,fullfile(fasSaveDir,sprintf('%s_cst_%s',mfilename,fasCleanedSaveName)))
close(h)

% Show the fasciles
h = figure;
% Get a random subset of the fibers
feConnectomeDisplay(arcuateFG,h,[.5 .7 .9]); hold on
feConnectomeDisplay(cstFG,    h,[.9 .7 .5]);
view(-80,0); camlight right; drawnow
feMakeMovie(h,numRotations,fullfile(fasSaveDir,sprintf('%s_arcuate_cst_%s',mfilename,fasCleanedSaveName)))
close(h)

end