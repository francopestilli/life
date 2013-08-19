function feFileToLoad = s_ms_test_connectomes_signal_prediction(connectomeType,rep)
%
% feFileToSave = s_ms_test_connectomes_signal_prediction(connectomeType,rep)
%
% Loads a connectome and produces a few scatter figures:
%  - the measured signal in data set 1 an data set 1 
%  - the LIFE model predicted signal and the measured signal.
%  - the Default* model predicted signal and the measured signal.
%
%  * Where the Default model is a model that assigns the same weight to all the fascicles. 
%
% Written by Franco Pestilli (c) Stanford University 2013 
if notDefined('connectomeType'), connectomeType = [4];end
if notDefined('rep'),            rep            = [1,2,3];end
if notDefined('cullType'),   cullType={'','culledL2'};end % We use the unculled connecotme
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_prediction');end
if notDefined('percentiles'), percentiles = {[0 100],[20 80],[10,90],[0 100],[0 50],[50 100],[0:5:100],[0 .25 .5 1 2 4 8 16 32 64 100]};end
if notDefined('makeScattePlot'), make2DmapPlot = true;end
if notDefined('doscatter'), doscatter = false;end

colors = {[.5 .5 .5],[.6 .6 .6],[.75 .75 .75]};
fontSiz = 16;

for icull = 1:length(cullType)
    for ip = 1:length(percentiles)
        x = [];y = [];  % We have three reps of each conenctome.
        for irep = 1:length(rep)
            [trackingType, lmax, bval, diffusionModelParams] = getConditions(connectomeType);
            [feFileToLoad, fName] = msBuildFeFileName(trackingType,lmax,bval,rep(irep),diffusionModelParams,cullType{icull});
            
            if (ip == 1)
                if exist(feFileToLoad,'file')
                    fprintf('[%s] Found connectome, loading it:\n%s\n',mfilename,feFileToLoad);
                    load(feFileToLoad)
                else
                    error('[%s] Cannot find the FE structure...\n%s\n',mfilename,feFileToLoad);
                end
            end
            
            for iPrctl = 1:length(percentiles{ip})-1
                drmse = feGetRep(fe,'vox rmse ratio');
                bins  = prctile(drmse,percentiles{ip});
                x(irep,:) = bins(1:end-1) + diff(bins)/2;
                
                % Find voxels where the Rrmse is less then 1
                goodvoxels = find( (drmse >= bins(iPrctl)) & (drmse < bins(iPrctl+1)));
                
                % Extract the data to dat reliability:
                LIFE_psig = feGet(fe, 'psigfiber',goodvoxels);
                Data1     = feGet(fe, 'dsigdemeaned',goodvoxels);
                Data2     = feGetRep(fe, 'dsigdemeaned',goodvoxels);
                
                % Compute the default model prediction
                DFM           = sum(feGet(fe,'Mfiber',goodvoxels),2);
                [~, DFMw, R2] = feFitModel(DFM,Data1,'lsqnonneg'); 
                DFM_psig      = full(DFM) * DFMw;
                
                if irep == 1
                    % LIFE Model prediction signal.
                    x_label='Data_{demeaned diffusion signal}';
                    y_label='LIFE model_{demeaned diffusion signal}';
                    if doscatter
                        figName = sprintf('LIFE_model_prediction_scatter%s_%s_prc_m%iM%i',fName,cullType{icull},100*percentiles{ip}(iPrctl),100*percentiles{ip}(iPrctl+1));
                        makeScatter(Data2,LIFE_psig,x_label, y_label, figName,colors{irep},fontSiz,fullfile(saveDir,'scatter'));
                    end
                    if make2DmapPlot
                        figName = sprintf('LIFE_model_prediction_2DMAP%s_%s_prc_m%iM%i',fName,cullType{icull},100*percentiles{ip}(iPrctl),100*percentiles{ip}(iPrctl+1));
                        [fh,sh] = make2Dmap(Data2,LIFE_psig,x_label, y_label, figName,fontSiz,fullfile(saveDir,'maps'),0);
                        figName = sprintf('LIFE_model_prediction_2DCON%s_%s_prc_m%iM%i',fName,cullType{icull},100*percentiles{ip}(iPrctl),100*percentiles{ip}(iPrctl+1));
                        [fh,sh] = make2Dmap(Data2,LIFE_psig,x_label, y_label, figName,fontSiz,fullfile(saveDir,'contours'),1);
                    end
                    
                    % Implicit Model prediction.
                    x_label='Data_{demenade diffusion signal}';
                    y_label='Default model_{demeaned diffusion signal}';
                    if doscatter
                        figName = sprintf('DFM_model_prediction_scatter%s_%s_prc_m%iM%i',fName,cullType{icull},100*percentiles{ip}(iPrctl),100*percentiles{ip}(iPrctl+1));
                        makeScatter(Data2,DFM_psig,x_label, y_label, figName,colors{irep},fontSiz,fullfile(saveDir,'scatters'));
                    end
                    if make2DmapPlot
                        figName = sprintf('DFM_model_prediction_2DMAP%s_%s_prc_m%iM%i',fName,cullType{icull},100*percentiles{ip}(iPrctl),100*percentiles{ip}(iPrctl+1));
                        [fh,sh] = make2Dmap(Data2,DFM_psig,x_label, y_label, figName,fontSiz,fullfile(saveDir,'maps'),0);
                        figName = sprintf('DFM_model_prediction_2DCON%s_%s_prc_m%iM%i',fName,cullType{icull},100*percentiles{ip}(iPrctl),100*percentiles{ip}(iPrctl+1));
                        [fh,sh] = make2Dmap(Data2,DFM_psig,x_label, y_label, figName,fontSiz,fullfile(saveDir,'contours'),1);
                    end
                    
                    % Data-to-data prediction.
                    x_label='Data_{one} (demeaned diffusion signal)';
                    y_label='Data_{two} (demeaned diffusion signal)';
                    if doscatter
                        figName = sprintf('Data2data_prediction_scatter%s_%s_prc_m%iM%i',fName,cullType{icull},100*percentiles{ip}(iPrctl),100*percentiles{ip}(iPrctl+1));
                        makeScatter(Data2,Data1,x_label, y_label, figName,colors{irep},fontSiz,fullfile(saveDir,'scatters'));
                    end
                    if make2DmapPlot
                        figName = sprintf('Data2data_prediction_2DMAP%s_%s_prc_m%iM%i',fName,cullType{icull},100*percentiles{ip}(iPrctl),100*percentiles{ip}(iPrctl+1));
                        [fh,sh] = make2Dmap(Data2,Data1,x_label, y_label, figName,fontSiz,fullfile(saveDir,'maps'),0);
                        figName = sprintf('Data2data_prediction_2DCON%s_%s_prc_m%iM%i',fName,cullType{icull},100*percentiles{ip}(iPrctl),100*percentiles{ip}(iPrctl+1));
                        [fh,sh] = make2Dmap(Data2,Data1,x_label, y_label, figName,fontSiz,fullfile(saveDir,'contours'),1);
                    end
                end
                
                % Store the correlation coefficient to make a nice plot across
                % Rrmse thresholds
                R.life(irep,iPrctl) = corr2(LIFE_psig,Data2);
                R.dfm(irep,iPrctl)  = corr2(DFM_psig,Data2);
                R.data(irep,iPrctl) = corr2(Data1,Data2);
                
                % Store the RMSE and RMSE ratio.
                rmse.m(irep,iPrctl) = median(feGetRep(fe,'vox rmse',      goodvoxels));
                rmse.d(irep,iPrctl) = median(feGetRep(fe,'vox rmse data', goodvoxels));
                rmse.r(irep,iPrctl) = median(feGetRep(fe,'vox rmse ratio',goodvoxels));
            end
        end
    % Make a plot of the increase in Correlation coefficient across Rrmse thresholds
    figName = sprintf('Prediction_improvement_Rrmse_%s_%s_prctGroup%i',fName,cullType{icull},ip);
    x_label = 'R_{rmse}';
    makePrediction(mean(rmse.r,1),R,x_label, y_label, figName,colors,fontSiz,saveDir);

    end
    
    
end
end

%-----------------------------------%
function [fh,sh] = make2Dmap(x,y,x_label, y_label, figName,fontSiz,saveDir,docontour)
%
% Make a 2D histogram (map), similar to the a scatter plot
%
% Franco Pestilli (c) 2013 Stanford University
fh = mrvNewGraphWin(figName);

% Makes a 2D density map.
ax.minmax = [-300 300];
ax.res   = 101;
ax.bins  = linspace(ax.minmax(1),ax.minmax(2),ax.res);
[ymap,xt] = hist3([x,y],{ax.bins, ax.bins});
ymap = (log10(ymap'));
xt{2} = (xt{2});
if ~docontour
    sh = imagesc(flipud(ymap));
    colormap(flipud(hot))
    hold on
    plot([0 ax.res],[ax.res 0],'k-');
else
    [~,sh] = contour(ymap,'levelstep',.25, 'linewidth',1);
    colormap(hot)
    hold on
    plot(ax.minmax,ax.minmax,'k-');
end
set(gca,'xlim',[0 ax.res],...
    'ylim',    [0 ax.res], ...
    'ytick',ceil([0 ax.res/2 ax.res]), ...
    'xtick',ceil([0 ax.res/2 ax.res]), ...
    'xticklabel',[xt{1}(1) xt{1}(ceil(end/2)) xt{1}(end)], ...
    'yticklabel',[xt{2}(1) xt{2}(ceil(end/2)) xt{2}(end)], ...
    'tickdir','out','box','off', ...
    'fontsize',fontSiz,'visible','on');
axis('square')
ylabel(y_label,'fontsize',fontSiz);
xlabel(x_label,'fontsize',fontSiz);
c = colorbar();
ytc = round(10.^get(c,'ytick'));
set(c,'yticklabel',ytc, ...
    'tickdir','out', ...
    'fontsize',fontSiz,'visible','on');
title(sprintf('Pearson correlation coefficient: %2.3f',corr2(x,y)));
saveFig(fh,fullfile(saveDir,figName),1);
end

%-----------------------------------%
function fh = makePrediction(x,R,x_label, y_label, figName,colors,fontSiz,saveDir)

% Sort the values
[x,idx] = sort(x);
life    = mean(R.life,1);life = life(idx);
data    = mean(R.data,1);data = data(idx);
dfm     = mean(R.dfm,1);dfm   = dfm(idx);

fh = mrvNewGraphWin(figName);
plot(x,life,'o-', ...
    'color',colors{1}, ...
    'markerfacecolor',colors{1}, ...
    'markeredgecolor','w', ...
    'markersize',12);
hold on
% Data reliability
plot(x,data,'o-', ...
    'color',colors{2}, ...
    'markerfacecolor',colors{2}, ...
    'markeredgecolor','w', ...
    'markersize',12);
% Default model
plot(x,dfm,'o-', ...
    'color',colors{3}, ...
    'markerfacecolor',colors{3}, ...
    'markeredgecolor','w', ...
    'markersize',12);
legend(gca,{'LIFE','Data','DFM'});
plot([x;x],[life;life] + [-std(R.life(:,idx),[],1);std(R.life(:,idx),[],1)],'r-', 'linewidth',2);
plot([x;x],[data;data] + [-std(R.data(:,idx),[],1);std(R.data(:,idx),[],1)],'r-', 'linewidth',2);
plot([x;x],[dfm;dfm]   + [-std(R.dfm(:,idx), [],1);std(R.dfm(:,idx), [],1)],'r-', 'linewidth',2);
set(gca,'xlim',[min(x)-.1 max(x)+.1],'ylim',[0.2 1], ...
    'ytick',[.2 .5 .8 1],'xtick',[.5 1 2], ...
    'tickdir','out','box','off', ...
    'fontsize',fontSiz,'xscale','log');
ylabel(y_label,'fontsize',fontSiz);
xlabel(x_label,'fontsize',fontSiz);
saveFig(fh,fullfile(saveDir,figName),1);
end
    
%-----------------------------------%
function fh = makeScatter(x,y,x_label, y_label, figName,colors,fontSiz,saveDir)
%
% Make a scatter plot of two diffusion signals
%
% Franco Pestilli (c) 2013 Stanford University 

fh = mrvNewGraphWin(figName);
plot([-300 300],[-300 300],'k-');
hold on
plot([median(x) median(y)],[-300 300],'k--');
plot([-300 300],[median(y) median(y)],'k--');
plot(x,y,'ro','color',colors,...
    'markerfacecolor',colors);
axis('square');
set(gca,'ylim',[-300 300],'xlim',[-300 300], ...
    'ytick',[-300 0 300],'xtick',[-300 0 300], ...
    'tickdir','out','box','off', ...
    'fontsize',fontSiz);
ylabel(y_label,'fontsize',fontSiz);
xlabel(x_label,'fontsize',fontSiz);
title(sprintf('Pearson correlation coefficient: %2.3f',corr2(x,y)));
saveFig(fh,fullfile(saveDir,figName),1);

end

%----------------------------------%
function saveFig(h,figName,eps)
%
% Saves a figure with publication quality
%
% Franco Pestilli (c) 2013 Stanford University
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end
    fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

if ~eps
    eval(sprintf('print(%s, ''-djpeg99'', ''-opengl'', ''%s'');', num2str(h),figName));
else
    eval(sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'');', num2str(h),figName));
end

end

%-------------------------------------%
function [algo, lmax, bval, diffModParamsType] = getConditions(runType)
%
% trackingType,lmax,bval,rep,diffusionModelParams
%
% Returns the desired options for a tractography.
% 
% Franco Pestilli (c) 2013 Stanford University

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
