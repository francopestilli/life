function feFileToSave = s_ms_cull_connectomes_figure_plot(connectomeType,rep)
%
% feFileToSave = s_ms_cull_connectomes_plot(connectomeType,rep)
%
% Loads a connectome and culls it down to the minimum number of fibers.
%
% Written by Franco Pestilli (c) Stanford University 2013 
if notDefined('rep'), rep = [1,2,3];end

if notDefined('fig_saveDir'), 
    fig_saveDir = fullfile('/home/frk/Dropbox','connectome_culling_figures');
    if ~exist(fig_saveDir,'dir');mkdir(fig_saveDir);end
end
fontSiz = 16;
numBins = 60;
%numIterToanalyze = 25;
iterToanalyze = [1 2 4 6 8 10 12 14 16 18 20 22 24];

% We have three reps of each conenctome.
for irep = 1:length(rep)
    [trackingType, lmax, bval, diffusionModelParams] = getConditions(connectomeType);
    
    % Information on the path to the files to load.
    % This is where the inputs will be loaded from
    feFileToLoad = msBuildFeFileName(trackingType,lmax,bval,rep(irep),diffusionModelParams);
    % This is the file where the output connectome will be saved.
    feFileToLoad = [feFileToLoad(1:end-4),'culledL2_example12thP','.mat'];
    
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
    percentReduxFibers(irep,:) = 100*(1 - numFibers(irep,:)./numFibers(irep,1));
end

% Root mean square error
figName = 'culling_cross_validated_rmse';
fh = mrvNewGraphWin(figName);
m_rmse  = mean(rmse,1);
x =1:length(m_rmse);
err_rmse = [m_rmse; m_rmse] + [-2*std(rmse); 2*std(rmse)];
plot([x;x],err_rmse,'r-','linewidth',3)
hold on
plot(x,m_rmse,'ko-', ...
    'color','k', ...
    'markerfacecolor','k',...
    'markeredgecolor','w', ...
    'markersize',8)
ylabel('cross-validated rmse')
xlabel('Iteration number')
%axis('square')
set(gca,'tickdir','out','box','off','FontSize',fontSiz,'xtick',iterToanalyze)
saveFig(fh,fullfile(fig_saveDir,figName))

% Root-mean squared ratio
figName = 'culling_cross_validated_Rrmse';
fh = mrvNewGraphWin(figName);
x =1:length(m_rmse);
m_rrmse  = mean(rrmse,1);
err_rrmse = [m_rrmse; m_rrmse] + [-2*std(rrmse); 2*std(rrmse)];
plot([x;x],err_rrmse,'r-','linewidth',3)
hold on
plot(x,m_rrmse,'ko-', ...
    'color','k', ...
    'markerfacecolor','k',...
    'markeredgecolor','w', ...
    'markersize',8)
ylabel('R_{rmse}')
xlabel('Iteration number')
%axis('square')
set(gca,'tickdir','out','box','off','FontSize',fontSiz,'xtick',iterToanalyze)
saveFig(fh,fullfile(fig_saveDir,figName))

% PErcent reduction in number of fibers from the original number
figName = 'culling_reduction_in_number_of_fibers';
fh = mrvNewGraphWin(figName);
x =1:length(percentReduxFibers);
m_prnf  = mean(percentReduxFibers,1);
err_prnf = [m_prnf; m_prnf] + [-2*std(percentReduxFibers); ....
                                 2*std(percentReduxFibers)];
plot([x;x],err_prnf,'r-','linewidth',3)
hold on
plot(x,m_prnf,'ko-', ...
    'color','k', ...
    'markerfacecolor','k',...
    'markeredgecolor','w', ...
    'markersize',8)
ylabel('Percent reduction in number of fascicles')
xlabel('Iteration number')
%axis('square')
set(gca,'tickdir','out','box','off','FontSize',fontSiz,'xtick',iterToanalyze)
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
ylabel('Number of fascicles')
xlabel('Iteration number')
set(gca,'tickdir','out','box','off','FontSize',fontSiz,'xtick',iterToanalyze)
saveFig(fh,fullfile(fig_saveDir,figName))

% Histograms of the weights across iterations:
maxWeight = 0.16;
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
colors = {[.65 .65 .65],[.35 .35 .35],[.35 .35 .35],[.35 .35 .35],[.35 .35 .35],[.35 .35 .35],[.35 .35 .35],[.35 .35 .35],...
          [.35 .35 .35],[.35 .35 .35],[.35 .35 .35],[.35 .35 .35],[.35 .35 .35],[.35 .35 .35],[.35 .35 .35],[.35 .35 .35], ...
          [.35 .35 .35],[.35 .35 .35],[.35 .35 .35],[.65 .65 .65],[.65 .65 .65],[.65 .65 .65], ...
          [.65 .65 .65],[.65 .65 .65],[.65 .65 .65]};
for iter= 1:length(iterToanalyze)
    figName = sprintf('weights_across_culling_iter%i',iterToanalyze(iter));
    fh = mrvNewGraphWin(figName);
    if iterToanalyze(iter) > 12 && iterToanalyze(iter) < 16
        ylims =  [0,0.14];
        yticks = [0,0.07,0.14];
    elseif iterToanalyze(iter) >= 16
            ylims =  [0,0.07];
        yticks = [0,0.035,0.07];
    else
        ylims =  [0,0.28];
        yticks = [0,0.07,0.14,0.21,.28];
    end
    
    title(sprintf('RMSE: %2.3f',cullingHist(1).rmsexv(iterToanalyze(iter))))
    set(gca,'tickdir','out','box','off','FontSize',fontSiz, ...
        'xlim',[0,maxWeight/2],'xtick',[0,maxWeight/4,maxWeight/2],...
        'ylim',ylims,'ytick',yticks,...
        'PlotBoxAspectRatio',[.5 1 1])
    hold on
    bar(xh,my{iter},'edgecolor',colors{iterToanalyze(iter)},'facecolor',colors{iterToanalyze(iter)})
    plot([xh;xh],[my{iter};my{iter}] + [-sdy{iter};sdy{iter}],'r-','linewidth',2)
    ylabel('Likelihood')
    xlabel('Fascicle contribution')
    saveFig(fh,fullfile(fig_saveDir,figName))
end
keyboard

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
function [algo, lmax, bval, diffModParamsType] = getConditions(runType)
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
