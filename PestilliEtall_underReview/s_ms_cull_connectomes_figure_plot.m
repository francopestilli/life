function feFileToSave = s_ms_cull_connectomes_figure_plot(connectomeType,rep)
%
% feFileToSave = s_ms_cull_connectomes_plot(connectomeType,rep)
%
% Loads a connectome and culls it down to the minimum number of fibers.
%
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

if notDefined('connectomeType'), connectomeType = [5];end
if notDefined('rep'), rep = [1,2];end
if notDefined('percentRemovedFas'),
    percentRemovedFas = 'culledL2_example15thP';
end

if notDefined('fig_saveDir'), 
    fig_saveDir = fullfile('/home/frk/Dropbox','connectome_culling_figures',percentRemovedFas);
    if ~exist(fig_saveDir,'dir');mkdir(fig_saveDir);end
end
fontSiz = 16;
numBins = 120;
iterToanalyze = [1 14];

% We have three reps of each conenctome.
for irep = 1:length(rep)
    [trackingType, lmax, bval, diffusionModelParams] = getConditionsLocal(connectomeType);
    
    % Information on the path to the files to load.
    % This is where the inputs will be loaded from
    feFileToLoad = msBuildFeFileName(trackingType,lmax,bval,rep(irep),diffusionModelParams);
    % This is the file where the output connectome will be saved.
    feFileToLoad = [feFileToLoad(1:end-4),percentRemovedFas,'.mat'];
    
    % Nowif this connectome was culled already we skip it.
    if exist(feFileToLoad,'file')
        fprintf('[%s] loading culled connectome:\n%s\n',mfilename,feFileToLoad);
        load(feFileToLoad,'cullingInfo')
    else
       error('[%s] cannot find culled connectome, skipping:\n%s\n',mfilename,feFileToLoad);
    end
    cullingHist(irep) = cullingInfo;
end

% Make aplot of rmse across iterations
for irep = 1:length(rep)
    rmse(irep,:)      = cullingHist(irep).rmsexv;
    rrmse(irep,:)     = cullingHist(irep).rrmse;
    numFibers(irep,:) = cullingHist(irep).numFibers;
    percentInitialNumberOfFascicles(irep,:) = 100.*(numFibers(irep,:)/numFibers(irep,1));
end

% Root mean square error
figName = 'culling_rmse_and_proportion_fibers';
fh = mrvNewGraphWin(figName);
% rmse 
rmse(:,1) = mean(rmse(:,2:5),2);
m_rmse    = mean(rmse,1);
err_rmse  = [m_rmse; m_rmse] + 2*[-std(rmse); std(rmse)];
x         = 1:length(m_rmse)-1;

% proportion of number of fibers
m_prnf  = mean(percentInitialNumberOfFascicles,1);

[ax,h1,h2] = plotyy(x,m_prnf(1:length(x)),x,m_rmse(1:length(x))); 
set(h2,'marker','s', ...
    'color','k', ...
    'markerfacecolor','k',...
    'markeredgecolor','w', ...
    'markersize',12)
set(h1,'Marker','o', ...
    'color',[.5 .5 .5], ...
    'markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','w', ...
    'markersize',10)
set(get(ax(2),'ylabel'),'String','Cross-validated rmse','FontSize',fontSiz);
set(ax(2), ...
    'tickdir','out','ycolor','k','box','off','FontSize',fontSiz, ...
    'xtick',[1 7 14 21 28],...
    'ytick',[26 28 30 32],'ylim',[25.875 32.125],'xlim',[0.5 29.5],...
        'PlotBoxAspectRatio',[2 1 1])
set(get(ax(1),'ylabel'),'String','Proportion fascicle number','FontSize',fontSiz);
set(ax(1),'ycolor',[.5 .5 .5],...
    'ylim',[-2 102],'tickdir','out','box','off','FontSize',fontSiz, ...
    'xtick',[1 7 14 21 28],...
    'ytick',[0 50 100],'xlim',[0.5 29.5],...
        'PlotBoxAspectRatio',[2 1 1])
axes(ax(2));hold on
plot([x;x],err_rmse(:,1:length(x)),'r-','linewidth',4)
xlabel('Iteration number','FontSize',fontSiz)
title(sprintf('Iter: [14 15] | %% fibers: [%i, %i] | rmse: [%2.3f, %2.3f] initial: %2.3f', ...
      round(m_prnf(14)),round(m_prnf(15)),m_rmse(14),m_rmse(15),m_rmse(1)))
saveFig(fh,fullfile(fig_saveDir,figName))

% Root-mean squared ratio
figName = 'culling_Rrmse_and_proportion_fibers';
fh = mrvNewGraphWin(figName);
x =1:length(m_rmse)-1;
rrmse(:,1) = mean(rrmse(:,2:5),2);
m_rrmse    = mean(rrmse,1);
err_rrmse  = [m_rrmse; m_rrmse] + [-2*std(rrmse); 2*std(rrmse)];
[ax,h1,h2] = plotyy(x,m_prnf(1:length(x)),x,m_rrmse(1:length(x))); 
set(h2,'marker','s', ...
    'color','k', ...
    'markerfacecolor','k',...
    'markeredgecolor','w', ...
    'markersize',12)
set(h1,'Marker','o', ...
    'color',[.5 .5 .5], ...
    'markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','w', ...
    'markersize',10)
set(get(ax(2),'ylabel'),'String','R_{rmse}','FontSize',fontSiz);
set(ax(2), ...
    'tickdir','out','ycolor','k','box','off','FontSize',fontSiz, ...
    'xtick',[1 7 14 21 28],...
    'ytick',[0.9 1 1.1],'ylim',[.9-0.005 1.1+0.005],'xlim',[0.5 29.5],...
        'PlotBoxAspectRatio',[2 1 1])
set(get(ax(1),'ylabel'),'String','Proportion fascicle number','FontSize',fontSiz);
set(ax(1),'ycolor',[.5 .5 .5],...
    'ylim',[-2 102],'tickdir','out','box','off','FontSize',fontSiz, ...
    'xtick',[1 7 14 21 28],...
    'ytick',[0 50 100],'xlim',[0.5 29.5],...
        'PlotBoxAspectRatio',[2 1 1])
axes(ax(2));hold on
plot([x;x],err_rrmse(:,1:length(x)),'r-','linewidth',4)
xlabel('Iteration number','FontSize',fontSiz);
title(sprintf('Iter: [14 15] | %% fibers: [%i, %i] | R_{rmse}: [%2.3f, %2.3f] initial: %2.3f', ...
      round(m_prnf(14)),round(m_prnf(15)),m_rrmse(14),m_rrmse(15),m_rrmse(1)))
saveFig(fh,fullfile(fig_saveDir,figName))

% Number of fibers across iteration
figName = 'culling_number_of_fibers';
fh = mrvNewGraphWin(figName);
x =1:length(numFibers);
m_nf  = mean(numFibers,1);
err_nf = [m_nf; m_nf] + [-2*std(numFibers); ....
                                 2*std(numFibers)];
plot([x;x],err_nf,'r-','linewidth',3)
hold on
plot(x,m_nf,'ko-', ...
    'color','k', ...
    'markerfacecolor','k',...
    'markeredgecolor','w', ...
    'markersize',8)
%axis('square')
ylabel('Number of fascicles','FontSize',fontSiz)
xlabel('Iteration number','FontSize',fontSiz)
set(gca,'tickdir','out','box','off','FontSize',fontSiz,'xtick',iterToanalyze,...
        'PlotBoxAspectRatio',[1 .5 1])
saveFig(fh,fullfile(fig_saveDir,figName))

% Histograms of the weights across iterations:
maxWeight = 0.08;
xh = linspace(0,maxWeight,numBins);
for iter= 1:length(iterToanalyze)
    for irep =1:length(rep)
        % We use always the same bins
        [yh{iter}(irep,:),~] = hist(cullingHist(irep).weights{iterToanalyze(iter)},xh);

        % Compute mean and std
        my{iter}  = mean(yh{iter}./sum(yh{iter}(:)));
        sdy{iter} = std(yh{iter}./sum(yh{iter}(:)));
    end
end

% Make a histogram of the weigths for eahc iteration
colors = {[.35 .35 .35]};
for iter= 1:length(iterToanalyze)
    figName = sprintf('weights_across_culling_iter%i',iterToanalyze(iter));
    fh = mrvNewGraphWin(figName);
    if iterToanalyze(iter) >= 8 && iterToanalyze(iter) < 12
        ylims =  [0,0.16];
        yticks = [0,0.08,0.16];
    elseif iterToanalyze(iter) >= 12
        ylims =  [0,0.04];
        yticks = [0,0.02,0.04];
    else
        ylims =  [0,0.4];
        yticks = [0,0.16,.32];
    end
    
    title(sprintf('RMSE: %2.3f',mean(rmse(:,iterToanalyze(iter)),1)))
    set(gca,'tickdir','out','box','off','FontSize',fontSiz, ...
        'xlim',[-0.005,maxWeight+0.005],'xtick',[0,maxWeight/2,maxWeight],...
        'ylim',ylims,'ytick',yticks,...
        'PlotBoxAspectRatio',[1 .5 1])
    hold on
    bar(xh,my{iter},'edgecolor',colors{1},'facecolor',colors{1})
    plot([xh;xh],[my{iter};my{iter}] + 2*[-sdy{iter};sdy{iter}],'r-','linewidth',2)
    ylabel('Probability','FontSize',fontSiz)
    xlabel('Fascicle weight','FontSize',fontSiz)
    saveFig(fh,fullfile(fig_saveDir,figName))
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)

printCommand = sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName);
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

% do the printing here:
eval(printCommand);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [algo, lmax, bval, diffModParamsType] = getConditionsLocal(runType)
%
% trackingType,lmax,bval,rep,diffusionModelParams
%
% Returns the desired options for a tractography.
switch runType
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 2000
    case 1
        algo = 'p';lmax = 2;bval=2000;diffModParamsType = [1,0];
    case 2
        algo = 'p';lmax = 4;bval=2000;diffModParamsType = [1,0];
    case 3
        algo = 'p';lmax = 6;bval=2000;diffModParamsType = [1,0];
    case 4
        algo = 'p';lmax = 8;bval=2000;diffModParamsType = [1,0];
    case 5
        algo = 'p';lmax = 10;bval=2000;diffModParamsType = [1,0];
    case 6
        algo = 'p';lmax = 12;bval=2000;diffModParamsType = [1,0];
    case 7
        algo = 'p';lmax = 14;bval=2000;diffModParamsType = [1,0];
    case 8
        algo = 'p';lmax = 16;bval=2000;diffModParamsType = [1,0];
      
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 4000
    case 9
        algo = 'p';lmax = 2;bval=4000;diffModParamsType = [1,0];
    case 10
        algo = 'p';lmax = 4;bval=4000;diffModParamsType = [1,0];
    case 11
        algo = 'p';lmax = 6;bval=4000;diffModParamsType = [1,0];
    case 12
        algo = 'p';lmax = 8;bval=4000;diffModParamsType = [1,0];
    case 13
        algo = 'p';lmax = 10;bval=4000;diffModParamsType = [1,0];
    case 14
        algo = 'p';lmax = 12;bval=4000;diffModParamsType = [1,0];
    case 15
        algo = 'p';lmax = 14;bval=4000;diffModParamsType = [1,0];
    case 16
        algo = 'p';lmax = 16;bval=4000;diffModParamsType = [1,0];
         
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 1000
    case 17
        algo = 'p';lmax = 2;bval=1000;diffModParamsType = [1,0];
    case 18
        algo = 'p';lmax = 4;bval=1000;diffModParamsType = [1,0];
    case 19
        algo = 'p';lmax = 6;bval=1000;diffModParamsType = [1,0];
    case 20
        algo = 'p';lmax = 8;bval=1000;diffModParamsType = [1,0];
    case 21
        algo = 'p';lmax = 10;bval=1000;diffModParamsType = [1,0];
    case 22
        algo = 'p';lmax = 12;bval=1000;diffModParamsType = [1,0];
    case 23
        algo = 'p';lmax = 14;bval=1000;diffModParamsType = [1,0];
    case 24
        algo = 'p';lmax = 16;bval=1000;diffModParamsType = [1,0];
         
    % MRTRIX deterministic tractogrpahy 
    % BVAL = 2000
    case 1+16
        algo = 'd';lmax = 2;bval=2000;diffModParamsType = [1,0];
    case 2+16
        algo = 'd';lmax = 4;bval=2000;diffModParamsType = [1,0];
    case 3+16
        algo = 'd';lmax = 6;bval=2000;diffModParamsType = [1,0];
    case 4+16
        algo = 'd';lmax = 8;bval=2000;diffModParamsType = [1,0];
    case 5+16
        algo = 'd';lmax = 10;bval=2000;diffModParamsType = [1,0];
    case 6+16
        algo = 'd';lmax = 12;bval=2000;diffModParamsType = [1,0];
    case 7+16
        algo = 'd';lmax = 14;bval=2000;diffModParamsType = [1,0];
    case 8+16
        algo = 'd';lmax = 16;bval=2000;diffModParamsType = [1,0];
      
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 4000
    case 9+16
        algo = 'd';lmax = 2;bval=4000;diffModParamsType = [1,0];
    case 10+16
        algo = 'd';lmax = 4;bval=4000;diffModParamsType = [1,0];
    case 11+16
        algo = 'd';lmax = 6;bval=4000;diffModParamsType = [1,0];
    case 12+16
        algo = 'd';lmax = 8;bval=4000;diffModParamsType = [1,0];
    case 13+16
        algo = 'd';lmax = 10;bval=4000;diffModParamsType = [1,0];
    case 14+16
        algo = 'd';lmax = 12;bval=4000;diffModParamsType = [1,0];
    case 15+16
        algo = 'd';lmax = 14;bval=4000;diffModParamsType = [1,0];
    case 16+16
        algo = 'd';lmax = 16;bval=4000;diffModParamsType = [1,0];
         
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 1000
    case 17+16
        algo = 'd';lmax = 2;bval=1000;diffModParamsType = [1,0];
    case 18+16
        algo = 'd';lmax = 4;bval=1000;diffModParamsType = [1,0];
    case 19+16
        algo = 'd';lmax = 6;bval=1000;diffModParamsType = [1,0];
    case 20+16
        algo = 'd';lmax = 8;bval=1000;diffModParamsType = [1,0];
    case 21+16
        algo = 'd';lmax = 10;bval=1000;diffModParamsType = [1,0];
    case 22+16
        algo = 'd';lmax = 12;bval=1000;diffModParamsType = [1,0];
    case 23+16
        algo = 'd';lmax = 14;bval=1000;diffModParamsType = [1,0];
    case 24+16
        algo = 'd';lmax = 16;bval=1000;diffModParamsType = [1,0];
        
    % MRTRIX tensor-based tractogrpahy 
    % BVAL = 2000,4000,1000
    case 25+16
        algo = 't';lmax = 2;bval=2000;diffModParamsType = [1,0];
    case 26+16
        algo = 't';lmax = 2;bval=4000;diffModParamsType = [1,0];
    case 27+16
        algo = 't';lmax = 2;bval=1000;diffModParamsType = [1,0];
  
  otherwise
        keyboard
end

end
