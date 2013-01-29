function s_ms_test_connectomes_plots(stats,trackingType,lmax,diffusionModelParams,recompute)
%
% Load FE structeres obtained by preprocessing connectomesconstrained within a
% region of interest and within the cortex and makes some basic plot of
% statistics in the connectomes
%
% See also:
%
% Franco (C) 2012 Stanford VISTA team.

if notDefined('trackingType'),trackingType = 'deterministic';end
if notDefined('lmax'),        lmax=[2 :2: 16];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1.8,0.2];end
if notDefined('recompute'), recompute=1;end
if notDefined('stats')
  % Type of statistcs to compute, we compute 2 each for the arcuate alone and
  % the arcuate plus the CST
  stats = {'vox rmse ratio','vox rmse','vox rmse ratio voxel wise','vox rmse voxel wise'};
end

% ROIs connectomes and saved paths
connectSubfolders = {'life_mrtrix_rep1','life_mrtrix_rep2','life_mrtrix_rep3'};
plots_saveDir     = fullfile('/home/frk/Dropbox','connectomes_plots');
projectDir        = '/azure/scr1/frk/150dirs_b1000_b2000_b4000';
loadDir           = fullfile(projectDir,'results');

% build a file name, for loading or saving the results.
lmaxs = num2str(lmax); lmaxs(isspace(lmaxs))='_';
fileName = fullfile(plots_saveDir,sprintf('%s_ax%srd%s_tracking%s_results_Feb2013_lmax%s.mat',mfilename,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)),upper(trackingType),lmaxs));

if recompute % Recompute the satistics
  % Initialize the results
  stat = nan(3,length(lmax),length(stats),length(connectSubfolders));
  
  for i_lmax = 1:length(lmax)
    % Get the connectome names
    switch trackingType
      case {'deterministic','d'}
        connectomeFile = { ...
          sprintf( '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax%i_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_stream-500000.pdb',lmax(i_lmax)), ...
          sprintf( '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax%i_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_stream-500000.pdb',lmax(i_lmax)),...
          sprintf( '0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_csd_lmax%i_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_brainmask_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_wm_stream-500000.pdb',lmax(i_lmax)),...
          };
      case {'probabilistic','p'}
        connectomeFile = { ...
          sprintf( '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax%i_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_prob-500000.pdb',lmax(i_lmax)), ...
          sprintf( '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax%i_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb',lmax(i_lmax)),...
          sprintf( '0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_csd_lmax%i_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_brainmask_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_wm_prob-500000.pdb',lmax(i_lmax));
          };
      otherwise
        keyboard
    end
    
    for i_bval = 1:length(connectomeFile)
      for irep = 1:length(connectSubfolders)
        
        % Set up file names
        % Find an identificative name for the connectome that is short enough:
        cName = [connectomeFile{i_bval}(1:57),'_',connectomeFile{i_bval}(end-16:end-4)];
        feLoadDir         = fullfile(loadDir,connectSubfolders{irep},'fe_structures');
        feLoadName        = sprintf('%s_diffModAx%sRd%s_%s',cName,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)));
        
        % File to load
        feFileToLoad =  fullfile(feLoadDir,[feLoadName,feLoadName,'.mat']);
        
        if (exist(feFileToLoad,'file') == 2)
          fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)
          load(feFileToLoad);
        else
          fprintf('[%s] FE file NOT found: \n%s\n ======================================== \n\n\n',mfilename,feFileToLoad)
          keyboard
        end
        
        
        % Compute some statistcs
        for is = 1:length(stats)
          stat(i_bval,i_lmax,is,irep)   = median(feGetRep(fe,  stats{is}));
        end
      end % irep
    end   % i_bval
  end     % i_lmax
  
  % Save the results collected so far
  save(fileName,'stat','stats')
else
  load(fileName,'stat','stats')
end

% Make plots
for ii = 1:length(stats)
  makeBarPlotAcrossLmax(stat, plots_saveDir,stats,ii,lmax,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)),trackingType)
end

end % Main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeBarPlotAcrossLmax(stat,plots_saveDir,stats,ii,lmax,ad,rd,trackingType)
% Make a plot across bval
colors = {[.5 .4 .65],[.6 .45 .4],[.45 .6 .35]};
xa     = {[1 1 1.5 1.5],[1.5 1.5 2 2],[2 2 2.5 2.5], [2.5 2.5 3 3], ...
          [3 3 3.5 3.5],[3.5 3.5 4 4],[4 4 4.5 4.5], [4.5 4.5 5 5]};
bval   =  [1000,2000,4000];

switch stats{ii}
  case {'vox rmse ratio'}
    ymin = 0.75;ymax = 1.1;
    ticks = [ymin 1 ymax];
    
  case {'vox rmse ratio voxel wise'}
    ymin = 0.75;ymax = 1.1;
    ticks = [ymin 1 ymax];
    
  case {'vox rmse'}
    ymin = 18;ymax = 36;
    ticks = [ymin (ymax+ymin)/2 ymax];
    
  case {'vox rmse voxel wise'}
    ymin = 18;ymax = 36;
    ticks = [ymin (ymax+ymin)/2 ymax];
    
  otherwise
    keyboard
end

for bb = 1:3
  figName = sprintf('%s_%s_ax%srd%s_bval%i_lmax_2to16',stats{ii}(~isspace(stats{ii})),trackingType,ad,rd,bval(bb));
  h = mrvNewGraphWin(figName);clf; hold on
  plot([0 7.5],[1 1],'k--')
  set(gca,'box','off','tickDir','out','ticklen',[0.014 0], ...
    'ylim',[ymin ymax],...
    'ytick',ticks, ...
    'xlim',[0 6],...
    'xtick',[], ...
    'xticklabel',{'Lmax=12','Lmax=2'})
  col = colors{bb};
  for i_lmax = 1:length(lmax)
    col = col*0.86;
    patch([xa{i_lmax}], [ymin mean(squeeze(stat(bb,i_lmax,ii,:))) mean(squeeze(stat(bb,i_lmax,ii,:))) ymin], ...
      'k','FaceColor',col,...
      'EdgeColor','w');
    plot([mean(xa{i_lmax}) mean(xa{i_lmax})],[mean(squeeze(stat(bb,i_lmax,ii,:)))-2*std(squeeze(stat(bb,i_lmax,ii,:))), ...
          mean(squeeze(stat(bb,i_lmax,ii,:)))+2*std(squeeze(stat(bb,i_lmax,ii,:)))], ...
      'r-','linewidth',8);
        ylabel(stats{ii})
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