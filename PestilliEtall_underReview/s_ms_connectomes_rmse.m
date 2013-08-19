function feFileToLoad = s_ms_connectomes_rmse(connectomeType,rep)
%
% feFileToSave = s_ms_connectomes_rmse(connectomeType,rep)
%
% Loads a connectome and produces a figure of the rmse M-D and rmse D-D.
%
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

if notDefined('connectomeType'), connectomeType = 6;end
if notDefined('rep'),         rep          = [1,2];end
if notDefined('cullType'),   cullType={'culledL2',''};end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_rmse');end
colors = {[.25 .25 .25],[.5 .5 .5],[.75 .75 .75]};
fontSiz = 16;
nBins= logspace(log10(.5),log10(2),25);%[.5:.05:2];

for icull = 1:length(cullType)
% We have three reps of each conenctome.
for irep = 1:length(rep)
    [trackingType, lmax, bval, diffusionModelParams] = getConditions(connectomeType);
    [feFileToLoad, fName] = msBuildFeFileName(trackingType,lmax,bval,rep(irep),diffusionModelParams,cullType{icull});
    
    if exist(feFileToLoad,'file')
        fprintf('[%s] Found connectome, loading it:\n%s\n',mfilename,feFileToLoad);
        load(feFileToLoad)
    else
        error('[%s] Cannot find the FE structure...\n%s\n',mfilename,feFileToLoad);
    end
    
    % Extract the data to dat reliability:
    Drmse = feGetRep(fe, 'vox rmse data');
    Mrmse = feGetRep(fe, 'vox rmse');
    Rrmse = feGetRep(fe, 'vox rmse ratio');
    
    % Make a scatter plot
    figNameRmse = sprintf('rmse_scatter_%s_%',fName,cullType{icull});
    if irep == 1,
        fhRmse = mrvNewGraphWin(figNameRmse);
        plot([0 90],[0 90],'k-')
        hold on
        plot([median(Drmse) median(Drmse)],[0 90],'k--')
        plot([0 90],[median(Mrmse) median(Mrmse)],'k--')
    end
    figure(fhRmse)
    plot(Drmse,Mrmse,'ro','color',colors{irep},...
        'markerfacecolor',colors{irep})
    axis('square')
    set(gca,'ylim',[0 90],'xlim',[0 90], ...
        'ytick',[0 45 90],'xtick',[0 45 90], ...
        'tickdir','out','box','off', ...
        'fontsize',fontSiz)
    ylabel('Model_{rmse}','fontsize',fontSiz)
    xlabel('Data_{rmse}','fontsize',fontSiz)
    saveFig(fhRmse,fullfile(saveDir,figNameRmse),1)
    
    if irep == 1
        % Make a 2D histogram pf the RMSEt
        figNameRmseMap = sprintf('rmse_map_hist_%s_%s',fName,cullType{icull});
        fhRmseMap = mrvNewGraphWin(figNameRmse);
        [ymap,x]  = hist3([Mrmse;Drmse]',{[1:.5:60], [1:.5:60]});
        ymap = ymap./length(Mrmse);
        sh = imagesc(flipud(log10(ymap)));
        cm = colormap(flipud(hot)); view(0,90);
        axis('square')
        set(gca,'xlim',[0 length(x{2})],...
            'ylim',    [0 length(x{2})], ...
            'ytick',[0 30 60 90 length(x{2})], ...
            'xtick',[0 30 60 90 length(x{2})], ...
            'yticklabel',[x{2}(end) 45 30 15 0], ...
            'xticklabel',[0 15 30 45 x{2}(end)], ...
            'tickdir','out','box','off', ...
            'fontsize',fontSiz,'visible','off')
        % Print a first version without the axis and drawins.
        saveFig(fhRmseMap,fullfile(saveDir,figNameRmseMap),0)
        
        % Print a second version with all the info.
        set(gca,'xlim',[0 length(x{2})],...
            'ylim',    [0 length(x{2})], ...
            'ytick',[0 30 60 90 length(x{2})], ...
            'xtick',[0 30 60 90 length(x{2})], ...
            'yticklabel',[x{2}(end) 45 30 15 0], ...
            'xticklabel',[0 15 30 45 x{2}(end)], ...
            'tickdir','out','ticklen',[.025 .05],'box','off', ...
            'fontsize',fontSiz','visible','on')
        hold on
        plot3([0 120],[120 0],[max(ymap(:)) max(ymap(:))],'k-','linewidth',1)
        ylabel('M_{rmse}','fontsize',fontSiz)
        xlabel('D_{rmse}','fontsize',fontSiz)
        cb = colorbar;
        tck = get(cb,'ytick');
        set(cb,'yTick',[min(tck)  mean(tck) max(tck)], ...
               'yTickLabel',round(10000*10.^[min(tck),...
                                 mean(tck), ...
                                 max(tck)])/10000, ...
               'tickdir','out','ticklen',[.025 .05],'box','on', ...
               'fontsize',fontSiz','visible','on')
        saveFig(fhRmseMap,fullfile(saveDir,figNameRmseMap),1);
    end

    % Compute probability of Rrmse'
    clear x 
    [y(irep,:),x] = hist(Rrmse,nBins);
    nSum = sum(y(irep,:));
    y(irep,:) = y(irep,:)./nSum;
end

% Plot the median and 2*std Rrmse across repeated tracking
my = median(y);
ey = [my;my] + [-std(y); std(y)];
   
% Make a histogram of the Rrmse
figNameRatio = sprintf('Rrmse_hist_%s_%s',fName,cullType{icull});
fhRatio = mrvNewGraphWin(figNameRatio);
%plot([median(Rrmse) median((Rrmse))],[0 .16],'k--')
px = x(x<=1);px =[px px(end)];
py = [my(x<=1) 0];
pp = patch(px,py,[.8 .8 .8],'edgecolor',[.8 .8 .8]);
hold on
plot([1 1],[0 .16],'k--')
plot(x,my,'o-','color',colors{1}, ...
    'markerfacecolor', colors{1}, ...
    'markeredgecolor','w',...
    'markersize',20,'linewidth',2)
plot([x;x],ey,'r-','linewidth',3)
set(gca,'tickdir','out','box','off', ...
    'fontsize',fontSiz,'ylim',[0 .16],'ytick',[0 .08 .16], ...
   'xtick',[.5 1 2],  'xscale','log')
ylabel('Probability','fontsize',fontSiz)
xlabel('R_{rmse}','fontsize',fontSiz)
title(sprintf('Proportion R_{rmse}<= 1: %2.3f',sum(my(x<=1))),'fontsize',fontSiz)
saveFig(fhRatio,fullfile(saveDir,figNameRatio),1)

end
end


function saveFig(h,figName,eps)
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end
    fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

if ~eps
    eval(sprintf('print(%s, ''-djpeg99'', ''-opengl'', ''%s'');', num2str(h),figName));
else
    eval(sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'');', num2str(h),figName));
end

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
