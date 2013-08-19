function s_ms_mt_ips(hemisphere,saveDir)
%
% This script performs a test of conenctivity of MT+ (LO1 and LO2) with
% IPS0. THe following are the steps we perform:
%  - It loads a culled half hemisphere connectome (FE structure). 
%  - It loads the MT+ ROI.
%  - It loads the IPS0 ROI. 
%  - It finds the fascicles connecting MT+ and IPS0
%  - It reduces the connectome to the voxels and fibers of the conenctions
%  - It Perform a bootstrap test WITH and WITHOUT the connection between
%    MT+ and IPS0.
%
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

% Handle parallel computing
if matlabpool('size') == 0
    c = parcluster;
    c.NumWorkers = 12;
    matlabpool(c);
end

if notDefined('saveDir'), 
    saveDir = sprintf('~/Dropbox/connectomes_plot_mt_ips/%s',mfilename);
end
if notDefined('hemisphere'), hemisphere = {'left','right'};end
if notDefined('plotSlices'), plotSlices = 1;end

nboots      = 10000;
nmontecarlo = 100;

for ih = 1:length(hemisphere)
% Set all the variables that depend on the hemisphere
switch hemisphere{ih}
    case {'left'}
        feFileName = 'fe_culled_FP_150_B2000_LMAX8_left.mat';
        ipsRoi     = 'LIPS0_cleaned.mat';
        mtRoi      = 'LTO1_O2_cleaned.mat';
        axisLims   = [-67 -18 -110 -40 -18 80];
        vw         = [-75,30];
        slices     = {[-18 0 0],[0 -40 0],[0 0 -14 ]};
        lght       = 'left';
        SLaxLims   = [-55 2 -120 120 -20 40 ];
        histcolor  = [0.4 0.4 0.4];

    case {'right'}
        feFileName = 'fe_culled_FP_150_B2000_LMAX8_right.mat';
        ipsRoi     = 'RIPS0_cleaned.mat';
        mtRoi      = 'RTO1_O2_cleaned.mat';
        axisLims   = [18 67 -110 -40 -18 80];
        vw         = [75,30];
        slices     = {[18 0 0],[0 -40 0],[0 0 -14 ]};
        lght       = 'right';
        SLaxLims   = [-2 55 -120 120 -20 40 ];
        histcolor  = [0 0 0];

    otherwise
        keyboard
end

% Load the FE structure
feDir      = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/fe_hemispheres';
fprintf('[%s] Loading the FE structure...\n',mfilename)
load(fullfile(feDir,feFileName));

% Load the ROIs
roiDir = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/mt_roi';
mtFileName   = fullfile(roiDir,mtRoi);
ips0FileName = fullfile(roiDir,ipsRoi);
mt   = dtiReadRoi(mtFileName);
ips0 = dtiReadRoi(ips0FileName);
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';

% Find the fascicles in the connectome that touch both MT+ and IPS0.
operations = {'and','and'};
rois       = {mt,ips0};
fg         = feGet(fe,'fg acpc');
[mtIpsFG, keepFascicles] = feSegmentFascicleFromConnectome(fg, rois, operations, 'mt_ips_zero');

% Clean the fibers by length, fibers that too long are likely to go far
% frontal and not just touch MT+ and IPS0.
[Lnorm, Lmm]   = mbaComputeFiberLengthDistribution(mtIpsFG, false);
maxSD          = 1; % Max standard deviation of the fibers to keep in the group.
fibers2delete  = Lnorm > maxSD;

% Now let's get the indices of the fibers in the FE structure:
fasIndices = find(keepFascicles);
fasIndices = fasIndices(fibers2delete);

% Now let's mark the fascicles as deleted.
keepFascicles(fasIndices) = false;

% Show te new fiber group
% High-resolution Anatomy
% t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';
% t1     = niftiRead(t1File);
% mtIpsFG = fgExtract(mtIpsFG,Lnorm < maxSD,'keep');
% figName = sprintf('Test_MT_IPS0_connection_brain_%s_%s',hemisphere{ih},feFileName(1:end-4));
% figureHandle = mrvNewGraphWin(figName);
% h  = mbaDisplayBrainSlice(t1, slices{1});
% hold on
% h  = mbaDisplayBrainSlice(t1, slices{2});
% h  = mbaDisplayBrainSlice(t1, slices{3});
% tractColor = [.8 .6 .2];%[.3 .7 .9];
% [figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(mtIpsFG.fibers,figureHandle,tractColor,'uniform');
% delete(lightHandle)
% view(vw(1),vw(2));
% axis(axisLims);
% lightHandle = camlight(lght);
% set(gcf,'Position',[0.0148 0.0148 .35 .87])
% drawnow
% saveFig(figureHandle,fullfile(saveDir,figName),'jpg')

% Reduce the connectome to the voxels and fibers of the connection between
% MT+ and IPS0.
[feWithoutFas, feWithFas, ~] = feTestFascicle(fe,keepFascicles,0);

if plotSlices
    slice  = {[0 -56 0],[0 -58 0],[0 -60 0],[0 -62 0],[0 -64 0],[0 -66 0], ...
        [0 -68 0],[0 -70 0],[0 -72 0],[0 -74 0]};
    t1     = niftiRead(t1File);
    
    for iS = 1:length(slice)
        % Make a figure of brain slice of the RMSE with the fascicle
        figName = sprintf('rmse_map_WITH_%s_slice%i_%s',feFileName(1:end-4),slice{iS}(find(slice{iS})),hemisphere{ih});
        shW = makeBrainMap(feWithFas,t1,slice{iS},SLaxLims,figName,saveDir);
        
        % Make a figure of brain slice of the RMSE without the fascicle
        figName = sprintf('rmse_map_WITHOUT_%s_slice%i_%s',feFileName(1:end-4),slice{iS}(find(slice{iS})),hemisphere{ih});
        shWO = makeBrainMap(feWithoutFas,t1,slice{iS},SLaxLims,figName,saveDir);
    end
end

% Make a plot of the R-squared
WITH.rmse      = mean(feGetRep(feWithFas,   'vox  rmse'));
WITH.rrmse     = median(feGetRep(feWithFas,   'vox  rmse ratio'));
WITH.rmseall   = (feGetRep(feWithFas,   'vox  rmse'));

WITHOUT.rmse    = median(feGetRep(feWithoutFas,'vox  rmse'));
WITHOUT.rrmse   = median(feGetRep(feWithoutFas,'vox  rmse ratio'));
WITHOUT.rmseall = (feGetRep(feWithoutFas,'vox  rmse'));

%% The following is the code for the bootstrap test on the MEAN rmse
e_rmse      = mean(WITHOUT.rmseall);
NullSet     = [WITH.rmseall];
sizeWith    = length(WITH.rmseall);

nullDistribution = nan(nboots,nmontecarlo);
for inm = 1:nmontecarlo
    parfor ibt = 1:nboots
        nullDistribution(ibt,inm) = mean(randsample(NullSet,sizeWith,true));
    end
    
    [y(:,inm),xhis] = hist(nullDistribution(:,inm),linspace(28.5,31.5,100));
    y(:,inm) = y(:,inm)./sum(y(:,inm));
end
y_m = mean(y,2);
y_e = [y_m, y_m] + 2*[-std(y,[],2),std(y,[],2)];

% Plot the null distribution and the empirical difference
figName = sprintf('Test_MT_IPS0_connection_rmse_mean_HIST_%s_%s',hemisphere{ih},feFileName(1:end-4));
fh = mrvNewGraphWin(figName);
%bar(x,y_m,'k')
%plot([xhis,xhis],y_e','color',[.3 .3 .3]); % Error bars for the distribution in
patch([min(xhis), xhis,xhis, max(xhis)],[0, y_e(:)', 0],histcolor,'FaceColor',histcolor,'EdgeColor',histcolor); % Distribution as the +/- 2SD
hold on
plot([e_rmse,e_rmse],[0 max(y_m)],'r-','linewidth',2)
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[0 0.125], ... 
        'xlim',[28,35], ...
        'ytick',[0 0.05 0.1], ...
        'xtick',[27 29 31 33 35], ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('rmse','fontsize',16')

% (3) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
if max(nullDistribution)<e_rmse
      p = 100*0.5/nboots;
else  p = sum(nullDistribution( nullDistribution > e_rmse ),1);
end
dprime_mean(:,ih) = repmat(e_rmse,1,nmontecarlo)./std(nullDistribution,[],1);
title(sprintf('The probability of obtaining the difference by chance is less than %2.6f%%\nStrength of connection evidence %2.3f',p,mean(dprime_mean(:,ih))), ...
    'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),'eps')

%% The following is the code for the bootstrap test on the MEDIAN rmse
e_rmse   = median(WITHOUT.rmseall);
NullSet  = [WITH.rmseall];
sizeWith = length(WITH.rmseall);

nullDistribution = nan(nboots,nmontecarlo);
for inm = 1:nmontecarlo
    parfor ibt = 1:nboots
        nullDistribution(ibt,inm) = median(randsample(NullSet,sizeWith,true));
    end
    
    [y(:,inm),xhis] = hist(nullDistribution(:,inm),linspace(28.5,31.5,100));
    y(:,inm) = y(:,inm)./sum(y(:,inm));
end
y_m = mean(y,2);
y_e = [y_m, y_m] + 2*[-std(y,[],2),std(y,[],2)];

% Plot the null distribution and the empirical difference
figName = sprintf('Test_MT_IPS0_connection_rmse_median_HIST_%s_%s',hemisphere{ih},feFileName(1:end-4));
fh = mrvNewGraphWin(figName);
%bar(x,y_m,'k')
%plot([xhis,xhis],y_e','color',[.3 .3 .3]); % Error bars for the distribution in
patch([min(xhis), xhis,xhis, max(xhis)],[0, y_e(:)', 0],histcolor,'FaceColor',histcolor,'EdgeColor',histcolor); % Distribution as the +/- 2SD
hold on
plot([e_rmse,e_rmse],[0 max(y_m)],'r-','linewidth',2)
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[0 0.125], ... 
        'xlim',[28,35], ...
        'ytick',[0 0.05 0.1], ...
        'xtick',[27 29 31 33 35], ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('rmse','fontsize',16')
% (3) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
if max(nullDistribution)<e_rmse
      p = 100*0.5/nboots;
else  p = sum(nullDistribution( nullDistribution > e_rmse ),1);
end
dprime_median(:,ih) = repmat(e_rmse,1,nmontecarlo)./std(nullDistribution,[],1);
title(sprintf('The probability of obtaining the difference by chance is less than %2.6f%%\nStrength of connection evidence %2.3f',p,mean(dprime_median(:,ih))), ...
    'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),'eps')


%% The following is the code for the bootstrap test of the difference in rmse
%
% Compute a test of the diference in rmse
% (1) Get the differece in rmse observed empiriclly
EmpiricalDiff = WITHOUT.rmse(1) - WITH.rmse(1);

% (2) Compute the Null distribution by:
% (2.1) Combine all the rmse from both WITH and WITHOUT.
% (2.2) Compute 10,000 distributions of rmse one for WITH one for WITHOUT
%       with the same numerosity of the initial samples
% (2.3) Compute the difference between the medians of these 10,000
%       distributions.
NullSet     = [WITHOUT.rmseall WITH.rmseall];
sizeWith    = length(WITH.rmseall);
sizeWithout = length(WITHOUT.rmseall);

clear y x y_m y_e ;

nullDistribution = nan(nboots,nmontecarlo);
for inm = 1:nmontecarlo
parfor ibt = 1:nboots
    % Y = RANDSAMPLE(POPULATION,K)
    BootWith    = median(randsample(NullSet,sizeWith,true));
    BootWithout = median(randsample(NullSet,sizeWithout,true));   
    nullDistribution(ibt,inm) = BootWithout - BootWith;
end
    
    [y(:,inm),xhis] = hist(nullDistribution(:,inm),linspace(-1.5,1.5,100));
    y(:,inm) = y(:,inm)./sum(y(:,inm));
end
y_m = mean(y,2);
y_e = [y_m, y_m] + 2*[-std(y,[],2),std(y,[],2)];

% Plot the null distribution and the empirical difference
figName = sprintf('Test_MT_IPS0_connection_rmseDiff_HIST_%s_%s',hemisphere{ih},feFileName(1:end-4));
fh = mrvNewGraphWin(figName);
%bar(x,y_m,'k')
%plot([xhis,xhis],y_e','color',[.3 .3 .3]); % Error bars for the distribution in
patch([min(xhis), xhis,xhis, max(xhis)],[0, y_e(:)', 0],histcolor,'FaceColor',histcolor,'EdgeColor',histcolor); % Distribution as the +/- 2SD
hold on
plot([EmpiricalDiff,EmpiricalDiff],[0 max(y_m)],'r-','linewidth',2)
axis([-2 2 0 0.1])
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[0 0.12], ...
        'xlim',[-2 3], ...
        'ytick',[0 0.05 0.1], ...
        'xtick',[-2 -1 0 1 2 3], ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('Difference in rmse','fontsize',16')

% (3) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
if max(nullDistribution)<EmpiricalDiff
      p = 100*1/nboots;
else  p = sum(nullDistribution( nullDistribution > EmpiricalDiff ));
end
dprime_d(:,ih) = repmat(EmpiricalDiff,1,nmontecarlo)./std(nullDistribution,[],1);
title(sprintf('The probability of obtaining the difference by chance is less than %2.6f%%\nStrength of connection evidence %2.3f',p,mean(dprime_d(:,ih))), ...
    'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),'eps')
 
% Make a scatter plot:
figName = sprintf('Test_MT_IPS0_connection_rmse_SCATTER_%s_%s',hemisphere{ih},feFileName(1:end-4));
fh = mrvNewGraphWin(figName);
hold on
plot([0 80],[0 80],'k-',[median(WITH.rmse),median(WITH.rmse)],[0 80],'k--',[0 80],[median(WITHOUT.rmse),median(WITHOUT.rmse)],'k--')
plot(WITH.rmseall,WITHOUT.rmseall,'o','color',[.35 .35 .95],'markerfacecolor',[.35 .35 .95],'MarkerSize',12);
axis square
set(gca,'tickdir','out','box','off', ...
        'xlim',[0 80], ...
        'ylim',[0 80], ...
        'ytick',[0 40 80], ...
        'xtick',[0 40 80], ...
        'ticklen',[0.025 0.05],'FontSize',16)
xlabel('rmse WITH Tract','fontsize',16)
ylabel('rmse WITHOUT Tract','fontsize',16)
saveFig(fh,fullfile(saveDir,figName),'eps')

figName = sprintf('Test_MT_IPS0_connection_rmse_MAP_%s_%s',hemisphere{ih},feFileName(1:end-4));
fh = mrvNewGraphWin(figName);
[ymap,x]= hist3([WITHOUT.rmseall;WITH.rmseall]',{[0:2:60], [0:2:60]});
ymap = ymap./length(WITHOUT.rmseall);
sh = imagesc(flipud(log10(ymap)));
cm = colormap(flipud(hot)); 
axis('square')
set(gca,'xlim',[0 length(x{2})+1],...
    'ylim',    [0 length(x{2})+1], ...
    'ytick',[0 (length(x{2})+1)/2 length(x{2})+1], ...
    'xtick',[0 (length(x{2})+1)/2 length(x{2})+1], ...
    'yticklabel',[x{2}(end) x{2}(round(end/2)) 0], ...
    'xticklabel',[0 x{2}(round(end/2)) x{2}(end)], ...
    'tickdir','out','box','off', ...
    'fontsize',16,'visible','off');
saveFig(fh,fullfile(saveDir,figName),'jpg')
hold on
plot([1 length(x{2})],[length(x{2}) 1],'k-');
cb = colorbar;
tck = get(cb,'ytick');
set(cb,'yTick',[min(tck)  mean(tck) max(tck)], ...
    'yTickLabel',round(1000*10.^[min(tck),...
    mean(tck), ...
    max(tck)])/1000, ...
    'tickdir','out','ticklen',[.025 .05],'box','on', ...
    'fontsize',16','visible','on')
set(gca,'visible','on')
saveFig(fh,fullfile(saveDir,figName),'eps')
end % Hemisphere

% Make a plot fo the strenght fo conenction evidnce:
figName = sprintf('Test_MT_IPS0_Sce_bar_%s_%s',hemisphere{ih},feFileName(1:end-4));
fh = mrvNewGraphWin(figName);
mStrength = mean(dprime_d,1);
eStrength = [mStrength;mStrength] + 2*[-std(dprime_d,[],1);std(dprime_d,[],1)];
hb = bar(mStrength,'facecolor','k');
hold on
plot([1 1; 2 2]',[eStrength],'r-','linewidth',4)
set(gca,'xlim',[0.5 2.5],...
    'ylim',    [0 7], ...
    'ytick',[0 3.5 7], ...
    'xtick',[1 2], ...
    'xticklabel',{hemisphere{1} hemisphere{2}}, ...
    'tickdir','out','box','off', ...
    'fontsize',16,'visible','on');
ylabel(sprintf('S_{ce}'),'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),'eps')

% Make a plot fo the strenght fo conenction evidnce:
figName = sprintf('Test_MT_IPS0_Sce_median_RMSE_bar_%s_%s',hemisphere{ih},feFileName(1:end-4));
fh = mrvNewGraphWin(figName);
mStrength = mean(dprime_median,1);
eStrength = [mStrength;mStrength] + 2*[-std(dprime_median,[],1);std(dprime_median,[],1)];
hb = bar(mStrength,'facecolor','k');
hold on
plot([1 1; 2 2]',[eStrength],'r-','linewidth',4)
set(gca,'xlim',[0.5 2.5],...
    'ylim',    [0 300], ...
    'ytick',[0 150 300], ...
    'xtick',[1 2], ...
    'xticklabel',{hemisphere{1} hemisphere{2}}, ...
    'tickdir','out','box','off', ...
    'fontsize',16,'visible','on');
ylabel(sprintf('S_{ce}'),'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),'eps')

% Make a plot fo the strenght fo conenction evidnce:
figName = sprintf('Test_MT_IPS0_Sce_mean_RMSE_bar_%s_%s',hemisphere{ih},feFileName(1:end-4));
fh = mrvNewGraphWin(figName);
mStrength = mean(dprime_mean,1);
eStrength = [mStrength;mStrength] + 2*[-std(dprime_mean,[],1);std(dprime_mean,[],1)];
hb = bar(mStrength,'facecolor','k');
hold on
plot([1 1; 2 2]',[eStrength],'r-','linewidth',4)
set(gca,'xlim',[0.5 2.5],...
    'ylim',    [0 300], ...
    'ytick',[0 150 300], ...
    'xtick',[1 2], ...
    'xticklabel',{hemisphere{1} hemisphere{2}}, ...
    'tickdir','out','box','off', ...
    'fontsize',16,'visible','on');
ylabel(sprintf('S_{ce}'),'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),'eps')

end % Main function

%%%%%%%%%%%%%%%%%%%%%%%
function [fh,sh] = makeBrainMap(fe,t1,slice,axLims,figName,saveDir)

% Make a map of the RMSE WITH and WITHOUT the fascicle:
coords  = feGet(fe,'roi coords') + 1;
xform   = feGet(fe,'xform img 2 acpc');
         
% Cross-validate RMSE
rmse = feGetRep(fe, 'vox rmse');
img  = feReplaceImageValues(nan(feGet(fe,'map size')),rmse,coords);
maxr = 50;

% Make anifti file from the rmse
ni  = niftiCreate('data',mbaNormalize(img,[0,1]), ...
    'qto_xyz',xform, ...
    'fname','rmse', ...
    'data_type',class(img));

% Open a figure
fh = mrvNewGraphWin(figName);

% Show the anatomy with the overlay
sh = mbaDisplayOverlay(t1, ni, slice, [], 'hot');

axis(axLims)

saveMap(fh,figName,saveDir,nanmean(img(:)),nanmedian(img(:)),nanstd(img(:)),maxr)
end

%---------------------------------%
function saveMap(fh,figName,saveDir,M,m,SD,maxfd)
% This helper function saves two figures for each map and eps with onlythe
% axis and a jpg with only the brain slice.
% The two can then be combined in illustrator.
%
% First we save only the slice as jpeg.
set(gca,'fontsize',16,'ztick',[-20 0 20 40], ...
    'xtick',[-50 -25 0 25 50], ...
    'tickdir','out','ticklength',[0.025 0])
axis off
saveFig(fh,fullfile(saveDir,'maps',figName),'tiff')
saveFig(fh,fullfile(saveDir,'maps',figName),'png')

% Then we save the slice with the axis as
% eps. This will only generate the axis
% that can be then combined in illustrator.
axis on
grid off

title(sprintf('mean %2.2f | median %2.2f | SD %2.2f', ...
    M,m,SD),'fontsize',16)
zlabel('Z (mm)','fontsize',16)
xlabel('X (mm)','fontsize',16)
cmap = colormap(hot(255));
colorbar('ytick',linspace(0,1,5),'yticklabel', ...    
    {linspace(0,1,5)*50}, ...
    'tickdir','out','ticklength',[0.025 0],'fontsize',16)
saveFig(fh,fullfile(saveDir,'maps',figName),1)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName,type)

% MAke sure the folder to save the figure exists
[p,f,e] = fileparts(figName);
[success,message] = mkdir(p);
if ~isempty(message), disp(sprintf('%s.',message));end

% Find out which type of figure and geenerate the proper printing command.
switch type
    case {0,'jpeg','jpg'}
        printCommand = (sprintf('print(%s, ''-djpeg90'',''-r500'' , ''-noui'', ''-opengl'', ''%s'')', num2str(h),figName));
    case {1,'eps'}
        printCommand = (sprintf('print(%s, ''-cmyk'', ''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName));
    case 'png'
        printCommand =  (sprintf('print(%s, ''-dpng'',''-r500'', ''%s'')', num2str(h),figName));
    case 'tiff'
        printCommand = (sprintf('print(%s, ''-dtiff'',''-r500'', ''%s'')', num2str(h),figName));
    case 'bmp'
        printCommand = (sprintf('print(%s, ''-dbmp256'',''-r500'', ''%s'')', num2str(h),figName));
    otherwise
        keyboard
end

% do the printing here:
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
eval(printCommand);
end