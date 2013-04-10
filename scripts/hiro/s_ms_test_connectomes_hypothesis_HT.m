function s_ms_test_connectomes_hypothesis_HT(feFileToLoad,saveDir,roisDir,roiNames,roiOperations)
%
% s_ms_test_connectomes_hypothesis(feFileToLoad,saveDir)
%
% Steps:
%   (1) Load one FE structure for the occipital connectome.
%   (2) Select a Fascicle and 
%   (3) Tests the importance of the fascicle within its white matter path
%
% This is part of a series of reproducible science scritps to be published
% with the LiFE mansucript.
%
% Written by Franco Pestilli (c) Stanford Vista Team 2013
if notDefined('feFileToLoad'), 
    error('[%s] Please pass afull-path to a culled FE structure...',mfilename);
end
if notDefined('saveDir'),
    error('[%s] Please pass a full-path to folder where to save the results...',mfilename);
end
if notDefined('roisDir'),
    error('[%s] Please pass a full-path to folder where the ROis are saved...',mfilename);
end
if notDefined('roiNames'),
    error('[%s] Please pass the names of the ROIs...',mfilename);
end
if notDefined('roiOperations'),
    error('[%s] Please pass the operation to apply for each roi passed in....',mfilename);
end

% These are the ROIs used to segment the fascicle
for iroi = 1:length(roiNames)
    rois{iroi} = fullfile(roisDir,roiNames);
end

% Get the fe structure
disp('loading the LiFE structure...')
if ischar(feFileToLoad)
    fprintf('Loading %s ...\n',feFileToLoad)
    load(feFileToLoad);
else
    fe  =feFileToLoad;
    clear feFileToLoad;
end

% Extract the fiber group from the FE structure
fg = feGet(fe,'fibers acpc');

% Perform a series of and and not operations to segement the Fascicle
[~, keepFascicles] = feSegmentFascicleFromConnectome(fg, rois, roiOperations, 'prob connectome');

% Remove all the voxels from the connectome except the ones where the
% fascicle passes through. Fit the new model.
[feWithoutFas, feWithFas, ~] = feTestFascicle(fe,keepFascicles,0);

% Compute the RMSE and Rrmse
WITH.rmse     = median(feGetRep(feWithFas,   'vox  rmse'));
WITH.rmseall  = (feGetRep(feWithFas,   'vox  rmse'));

WITHOUT.rmse     = median(feGetRep(feWithoutFas,'vox  rmse'));
WITHOUT.rmseall  = (feGetRep(feWithoutFas,'vox  rmse'));

% Compute a statistic test of the difference in rmse between the case with
% the fascicle and the case without the fascicle.
%
% (1) Get the differece in rmse observed empiriclly
EmpiricalDiff = WITHOUT.rmse - WITH.rmse;

% (2) Compute the Null distribution by:
% (2.1) Combine all the rmse from both WITH and WITHOUT.
% (2.2) Compute 1000 samples of rmse one for WITH one for WITHOUT
%       with the same numerosity of the initial samples
% (2.3) Compute the difference between the medians of these 10,000
%       distributions.
NullSet     = [WITHOUT.rmseall WITH.rmseall];
sizeWith    = length(WITH.rmseall);
sizeWithout = length(WITHOUT.rmseall);

nboots = 1000;
nullDistribution = nan(nboots,1);
parfor ibt = 1:nboots
    % Y = RANDSAMPLE(POPULATION,K)
    BootWith    = median(randsample(NullSet,sizeWith));
    BootWithout = median(randsample(NullSet,sizeWithout));   
    nullDistribution(ibt) = BootWithout - BootWith;
end

% (3) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
if max(nullDistribution)<EmpiricalDiff
    p = 100*1/nboots;
else
    p = sum(nullDistribution(sort(nullDistribution)>EmpiricalDiff));
end

% Plot the null distribution and the empirical difference
figName = sprintf('Test_V4_V3ab_rmse_Probability');
fh = mrvNewGraphWin(figName);
[y,x] = hist(nullDistribution,100);
y = y./sum(y);
bar(x,y,'k')
hold on
plot([EmpiricalDiff,EmpiricalDiff],[0 max(y)],'r-','linewidth',2)
set(gca,'tickdir','out','box','off', 'ylim',[0 max(y)],'FontSize',16)
ylabel('Probability')
xlabel('Difference in rmse')
title(sprintf('The probability of obtaining the difference by chance is less then %2.3f%%',p), ...
    'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),1)
 
% Make a plot of the R-squared
figName = sprintf('Test_V4_V3ab_Rrmse');
fh = mrvNewGraphWin('Test ILF removal');
bar([WITH.rrmsem,WITHOUT.rrmsem],'FaceColor',colors{1})
hold on
plot([1 1; 2 2]',[WITH.rrmsesd,WITHOUT.rrmsesd],'r-','linewidth',16)
plot([0 3],[1 1],'k-')
ylabel('R_{rmse}')
set(gca,'xticklabel',{'With fascicle','Without fascicle'},'tickdir','out','box','off', ...
    'ylim',[0.75 1.4],'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),1)

% Make a scatter plto:
figName = sprintf('Test_V4_V3ab_rmse_SCATTER');
fh = mrvNewGraphWin(figName);
hold on
plot([0 80],[0 80],'k-',[median(WITHOUT.rmse),median(WITHOUT.rmse)],[0 80],'k--',[0 80],[median(WITH.rmse),median(WITH.rmse)],'k--')
plot(WITHOUT.rmseall,WITH.rmseall,'o','color',colors{ii},'markerfacecolor',colors{ii});
axis square
set(gca,'tickdir','out','box','off','xlim',[0 80],'ylim',[0 80],'FontSize',16)
ylabel('rmse WITH Fascicle')
xlabel('rmse WITHOUT Fascicle')
saveFig(fh,fullfile(saveDir,figName),1)

end % END function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName,eps)
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end
    fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

if ~eps
    eval(sprintf('print(%s, ''-djpeg90'', ''-opengl'', ''%s'')', num2str(h),figName));
else
    eval(sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName));
end

end
