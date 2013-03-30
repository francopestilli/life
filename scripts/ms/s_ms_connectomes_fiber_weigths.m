function s_ms_connectomes_fiber_weigths(stats,trackingType,lmax,diffusionModelParams,recompute)
%
% Load FE structeres obtained by preprocessing connectomesconstrained within a
% region of interest and within the cortex and makes some basic plot of
% statistics in the connectomes
%
% See also:
%
% Franco (C) 2012 Stanford VISTA team.

if notDefined('trackingType'),trackingType = 'deterministic';end
if notDefined('lmax'),        lmax=[2:2:16];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('recompute'), recompute=1;end
if notDefined('stats')
  % Type of statistcs to compute, we compute 2 each for the arcuate alone and
  % the arcuate plus the CST
  stats = {'weight global reliability', 'weight length corr', ...
           ...'weigths along fiber in voxel wise fit', ...
           ...'weigths correlation mean voxel wise vs. global',...'weight voxel wise reliability', ...
           };
  % {'weight voxel wise reliability'};
  % {'weight global reliability'};
  % {'weight length corr'};
  % {'num zero weigths per fiber in voxel wise fit'};
  % {'weight voxel wise reliability'};
  % {'num zero weigths per fiber in voxel wise fit'};
  % {'weigths along fiber in voxel wise fit'};
  % {'weigths difference along fiber in voxel wise fit'}
end

% ROIs connectomes and saved paths
connectSubfolders = {'life_mrtrix_rep3'};%,'life_mrtrix_rep2','life_mrtrix_rep3'};
plots_saveDir     = fullfile('/home/frk/Dropbox','connectomes_plots');
projectDir        = '/azure/scr1/frk/150dirs_b1000_b2000_b4000';
loadDir           = fullfile(projectDir,'results');

% build a file name, for loading or saving the results.
lmaxs = num2str(lmax); lmaxs(isspace(lmaxs))='_';
fileName = fullfile(plots_saveDir,sprintf('%s_ax%srd%s_tracking%s_results_Feb2013_lmax%s.mat',mfilename,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)),upper(trackingType),lmaxs));

if recompute % Recompute the satistics
  % Initialize the results
  stat = nan(length(stats),3,length(lmax),length(connectSubfolders));
  
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
    
    for irep = 1:length(connectSubfolders)
      for i_bval = 1:length(connectomeFile)
        
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
        
        for is = 1:length(stats)
          % Compute some statistcs
          switch stats{is}
            case {'weight length corr'}
              disp('Generating plot: ''weight length corr''')
              % Correlation plot between fiber length and weight magnitude in
              % a Global Fit
              w = (feGet(fe,  'fiber weights'));
              l   = fefgGet(feGet(fe,'fg acpc'),'length');
              makeWeigthLengthCorrelationPlot(w,l, plots_saveDir,lmax,i_lmax,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)),trackingType,irep,i_bval);
              
            case {'weight global reliability'}
              % Plot the reliability of the weights in each fiberfor a global
              % fit (one weight per fiber)
              disp('Generating plot: ''weight global reliability''')
              nBoots = 10;
              w = nan(feGet(fe,'n fibers'),nBoots);
              R2 = nan(1,nBoots);
              for iboot = 1:nBoots
                [~, w(:,iboot), R2(iboot)] = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsig demeaned'),'sgdnn');
              end
              
              % Make a plot
              makeWeightGlobalReliabilityPlot(w, R2, plots_saveDir,lmax,i_lmax,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)),trackingType,irep,i_bval)
              
            case {'weight voxel wise reliability'}
              % Plot the reliability of the weights in each fiberfor a global
              % fit (one weight per fiber)
              disp('Generating plot: ''weight voxel wise reliability''')
              nBoots = 5;
              siz = size(feGet(fe,'fiber weights voxel wise'));
              mw = nan(siz(2),nBoots);
              for iboot = 1:nBoots
                [~, w] = feFitModelByVoxel(fe);
                mw(:,iboot) = mean(w.weights,1);
              end
                
              % Make a plot
              makeWeightVWReliabilityPlot(mw, plots_saveDir,lmax,i_lmax,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)),trackingType,irep,i_bval)
              
            case {'weigths along fiber in voxel wise fit'}
              disp('Generating plot: ''weigths along fiber in voxel wise fit''')
              if notDefined('wvw'),  wvw        = feGet(fe,'fiber weights voxel wise');     end
              if notDefined('wvw1'),[wvw1 mwvw] = computFiberWeightsAlongTrajectory(wvw,fe);end  
              plotFiberWeigthAlongTrajectory(wvw,wvw1,plots_saveDir,lmax,i_lmax,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)),trackingType,irep,i_bval)
              
            case {'weigths correlation mean voxel wise vs. global'}
              disp('Generating plot: ''weigths correlation mean voxel wise vs. global''')
              if notDefined('wvw'),  wvw = feGet(fe,'fiber weights voxel wise');end
              if notDefined('wg'),   wg  = feGet(fe,'fiber weights');end
              if notDefined('mwvw'),[wvw1 mwvw] = computFiberWeightsAlongTrajectory(wvw,fe);end
              makeWeightGlobalWoxelWiseCorrPlot(wg,mwvw,plots_saveDir,lmax,i_lmax,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)),trackingType,irep,i_bval)
                            
            otherwise
              keyboard
          end
        end
        clear mwvw wvw wg w wvw1
        
      end   % i_bval
    end     % i_lmax
  end       % irep
  
  % Save the results collected so far
  %save(fileName,'stat','stats')
else
  load(fileName,'stat','stats')
end
  
end % Main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotFiberWeigthAlongTrajectory(w,wvw1,plots_saveDir,lmax,i_lmax,ad,rd,trackingType,irep,bb)
nFibers = 4;      % number of fibers to plot the weights for
minW    = 0.0001; % minimum weight to plot.

% Default plot format
colors = {[.5 .4 .65],[.6 .45 .4],[.45 .6 .35]};
bval   =  [1000,2000,4000];
col = colors{bb};
col = col*0.86;
sym = {'ro-','rs-','rd-','r^-'};

% Find nFibers with large weights
wmean        = mean(w); % The mean here will be biased by large weights.
wMeanNotZero = find(wmean > minW);
w_index      = wMeanNotZero( randi(length(wMeanNotZero),nFibers,1) );

figName = sprintf('weights_along_fascile_voxel_wise_%s_ax%srd%s_bval%i_lmax%i_rep%i',trackingType,ad,rd,bval(bb),lmax(i_lmax),irep);
h = mrvNewGraphWin(figName);clf; hold on

for i_fib = 1:nFibers
  plot(wvw1{w_index(i_fib)},sym{i_fib},'color',col,'MarkerFaceColor',col)
  hold on
  xlabel('Voxel in fascicle path')
  ylabel('Fascicle contribution (Voxel-wise fit)')
end
set(gca,'xlim',[0 30], 'box','off','tickDir','out')

% Save the figure
saveFig(h,fullfile(plots_saveDir,figName))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wvw1 mwvw] = computFiberWeightsAlongTrajectory(wvw,fe)

% Get the number of fiers and the number of voxels
nFibers = size(wvw,2);
nVox = length(fe.life.voxel2FNpair);

wvw1 = cell(1,nFibers);
mwvw =  nan(1,nFibers);
swvw = mwvw;%mdwvw = swvw;
vx_index = false(1,nVox);

tic
for i_fib = 1:nFibers
  for ii = 1:nVox
    vx_index(ii) = any( fe.life.voxel2FNpair{ii}(:,1) == i_fib );
  end
  
  wvw1{i_fib}  =    wvw(vx_index,i_fib);
  mwvw(i_fib)  =   mean(wvw1{i_fib});
  %mdwvw(i_fib) = median(wvw1{i_fib});
  %swvw(i_fib)  =    sum(wvw1{i_fib});
end
toc

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeWeightGlobalWoxelWiseCorrPlot(wg,mwvw,plots_saveDir,lmax,i_lmax,ad,rd,trackingType,irep,bb)

colors = {[.5 .4 .65],[.6 .45 .4],[.45 .6 .35]};
bval   =  [1000,2000,4000];
col = colors{bb};
col = col*0.86;

figName = sprintf('weights_global_woxel_wise_corr_%s_ax%srd%s_bval%i_lmax%i_rep%i',trackingType,ad,rd,bval(bb),lmax(i_lmax),irep);
h = mrvNewGraphWin(figName);clf; hold on
thiscor = corr(wg(:),mwvw(:));
plot(wg,mwvw,'o','color',col)
axis square; set(gca,'box','off')
title(sprintf('Pearson correlation coefficient: %2.3f',thiscor))
xlabel('Fascicle contribution (global fit)')
ylabel('Fascicle mean contribution (voxel-wise fit)')
set(gca,'box','off','tickDir','out')

% Save the figure
saveFig(h,fullfile(plots_saveDir,figName))
end
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeWeightVWReliabilityPlot(w,plots_saveDir,lmax,i_lmax,ad,rd,trackingType,irep,bb)
% Make a plot across bval
colors = {[.5 .4 .65],[.6 .45 .4],[.45 .6 .35]};
bval   =  [1000,2000,4000];
nFibers = 200; % Number of fibers to test

figName = sprintf('weights_VW_reliability_%s_ax%srd%s_bval%i_lmax%i_rep%i',trackingType,ad,rd,bval(bb),lmax(i_lmax),irep);
h = mrvNewGraphWin(figName);clf; hold on
col = colors{bb};
col = col*0.86;

% Select a random set of fibers
w_indx = randi(size(w,1),nFibers,1);
% Plot the median weight values
[val, indx] = sort(median(w(w_indx,:),2),'descend');
std_val = std(w(w_indx,:),[],2);
ste_val = std_val(indx)/sqrt(length(indx));
hb = bar(val,1,'k');
hold on
% Plot error bars
x = get(hb,'XData');
plot([x;x],[val - ste_val, val + ste_val]', 'r-','LineWidth',3)
axis tight
ylabel('Fascicle contribution (a.u.)')
xlabel('Fascicle index')
set(gca,'box','off','tickDir','out')
drawnow

% Save the figure
saveFig(h,fullfile(plots_saveDir,figName))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeWeightGlobalReliabilityPlot(w,R2, plots_saveDir,lmax,i_lmax,ad,rd,trackingType,irep,bb)
% Make a plot across bval
colors = {[.5 .4 .65],[.6 .45 .4],[.45 .6 .35]};
bval   =  [1000,2000,4000];
nFibers = 200; % Number of fibers to test

figName = sprintf('weights_global_reliability_%s_ax%srd%s_bval%i_lmax%i_rep%i',trackingType,ad,rd,bval(bb),lmax(i_lmax),irep);
h = mrvNewGraphWin(figName);clf; hold on
col = colors{bb};
col = col*0.86;

% Select a random set of fibers
w_indx = randi(size(w,1),nFibers,1);
% Plot the median weight values
[val, indx] = sort(median(w(w_indx,:),2),'descend');
std_val = std(w(w_indx,:),[],2);
ste_val = std_val(indx)/sqrt(length(indx));
hb = bar(val,1,'k');
hold on
% Plot error bars
x = get(hb,'XData');
plot([x;x],[val - ste_val, val + ste_val]', 'r-','LineWidth',3)
axis tight
title(sprintf('R^2: mean = %2.3f, ste = %2.3f',mean(R2),std(R2)/sqrt(length(R2))))
ylabel('Fascicle contribution (a.u.)')
xlabel('Fascicle index')
set(gca,'box','off','tickDir','out')
drawnow

% Save the figure
saveFig(h,fullfile(plots_saveDir,figName))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeWeigthLengthCorrelationPlot(W,len,plots_saveDir,lmax,i_lmax,ad,rd,trackingType,irep,bb)
% Make a plot across bval
colors = {[.5 .4 .65],[.6 .45 .4],[.45 .6 .35]};
bval   =  [1000,2000,4000];

figName = sprintf('weights_length_corr_%s_ax%srd%s_bval%i_lmax%i_rep%i',trackingType,ad,rd,bval(bb),lmax(i_lmax),irep);
h = mrvNewGraphWin(figName);clf; hold on
col = colors{bb};
col = col*0.86;
thiscor = corr(W(:),len(:));
% Make plot of correlatio between length of fibers and weight of
% the fiber for the global weights
plot(W,len,'ok','color',col,'MarkerFaceColor',col,'MarkerSize',8)
axis square
set(gca,'xLim',[-0.001 max(W)+0.001],'TickDir','out','box','off')
title(sprintf('Pearson correlation coefficient: %2.3f',thiscor))
ylabel('Fascicle length (mm)')
xlabel('Fascicle contribution')
drawnow
% Save the figure
saveFig(h,fullfile(plots_saveDir,figName))
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