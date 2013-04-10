function s_ms_test_connectomes_hypothesis(feFileToLoad,trackingType,lmax,bval,rep)
%
% s_ms_test_connectomes_hypothesis(trackingType,lmax,bval,rep)
%
% Load a series of FE structures for conectomes of the occipital lobe.
% Removes the ILF/OR and tests the importance of the fascicle
%
% This is part of a series of reproducible science scritps to be published
% with the LiFE mansucript.
%
% Franco (c) Stanford Vista Team 2013
if notDefined('trackingType'),trackingType = 'p';end
if notDefined('lmax'),        lmax         = [2:2:16];end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = [1,2,3];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('cullType'),   cullType='culledL2';end
if notDefined('displayFascicles'), displayFascicles = 1;end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_ILF_subtration');end

% High-resolution Anatomy
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';

% Path to the ROIs used to segment the fascicle
roisPath = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/ROIs_restric_connectomes';

% These are the ROIs used to segment the fascicle
rois = {fullfile(roisPath,'ROI_Anterior_new.mat'),  ...
    fullfile(roisPath,'ROI_Inferior_new.mat'),  ...
    fullfile(roisPath,'ROI_posterior_new.mat'),  ...
    fullfile(roisPath,'ROI_Superior_new.mat')};

% These are the logical operations to apply to each ROI
operation = {'and','not','and','not'};

% Make some plots
colors = {[.5 .4 .65],[.6 .45 .4],[.45 .6 .35]};

for i_lmax = 1:length(lmax)
% Check that the connectome were proprocessed before attempting to make a
% plot.
done = s_ms_check_processes([],trackingType,lmax(i_lmax),bval,cullType);

if ~all(done), 
    fprintf('\n[%s] Not all connectomes were proprocessed... not proceeding with plot.\n',mfilename);
return
end

for irep = 1:length(rep)
        % Information on the path to the files to load.
        % This is where the inputs will be loaded from
        [feFileToLoad, fname] = msBuildFeFileName(trackingType,lmax(i_lmax),bval,rep(irep),diffusionModelParams,cullType);

    
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
    
    % Perform a series of end and not operations to segement the OR/ILF
    [fas, keepFascicles] = feSegmentFascicleFromConnectome(fg, rois, operation, 'prob connectome');
    
    if displayFascicles && (irep == 1)
        weights = feGet(fe,'fiber weights');
        figName = sprintf('Test_ILF_SUBTRACTION_Connectome_rep%i_%s_lmax%i',irep,fname,lmax(i_lmax));
        [~, lightHandle, sHcc] = mbaDisplayConnectome(fg.fibers,figure,[.8 .2 .1] ,'single',[],0.005+weights);
        delete(lightHandle)
        [~, lightHandle, sHfas] = mbaDisplayConnectome(fas.fibers,gcf,[] ,'single',[],.7);
        delete(lightHandle)
        lightHandle = camlight('right');
        set(gcf,'color','w')
        drawnow
        saveFig(gcf,fullfile(saveDir,figName),0)
        close all
    end
    
    % Remove all the voxels from the connectome except the ones where the
    % fascicle passes through. Fit the new model.
    [feWithoutFas, feWithFas, ~] = feTestFascicle(fe,keepFascicles,0);
    
    % Make a plot of the R-squared
    WITH.r2(irep)       = median(feGetRep(feWithFas,   'vox  r2'));
    WITH.rmse(irep)     = median(feGetRep(feWithFas,   'vox  rmse'));
    WITH.rrmse(irep)    = median(feGetRep(feWithFas,   'vox  rmse ratio'));
    WITH.rmseall{irep}  = (feGetRep(feWithFas,   'vox  rmse'));

    WITHOUT.r2(irep)    = median(feGetRep(feWithoutFas,'vox  r2'));
    WITHOUT.rmse(irep)  = median(feGetRep(feWithoutFas,'vox  rmse'));
    WITHOUT.rrmse(irep) = median(feGetRep(feWithoutFas,'vox  rmse ratio'));
    WITHOUT.rmseall{irep}  = (feGetRep(feWithoutFas,'vox  rmse'));

end

% Now make averages and std or the rmse, Rrmse, r2 values for plotting
WITH.r2m     = mean(WITH.r2);
WITH.r2sd    = [WITH.r2m-std(WITH.r2); WITH.r2m+std(WITH.r2)];
WITH.rmsem   = mean(WITH.rmse);
WITH.rmsesd  = [WITH.rmsem-std(WITH.rmse); WITH.rmsem+std(WITH.rmse)];
WITH.rrmsem  = mean(WITH.rrmse);
WITH.rrmsesd = [WITH.rrmsem-std(WITH.rrmse); WITH.rrmsem+std(WITH.rrmse)];

WITHOUT.r2m     = mean(WITHOUT.r2);
WITHOUT.r2sd    = [WITHOUT.r2m-std(WITHOUT.r2); WITHOUT.r2m+std(WITHOUT.r2)];
WITHOUT.rmsem   = mean(WITHOUT.rmse);
WITHOUT.rmsesd  = [WITHOUT.rmsem-std(WITHOUT.rmse); WITHOUT.rmsem+std(WITHOUT.rmse)];
WITHOUT.rrmsem  = mean(WITHOUT.rrmse);
WITHOUT.rrmsesd = [WITHOUT.rrmsem-std(WITHOUT.rrmse); WITHOUT.rrmsem+std(WITHOUT.rrmse)];

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
figName = sprintf('Test_ILF_SUBTRACTION_BOOT_test_rmse_%s_lmax%i',fname,lmax(i_lmax));
fh = mrvNewGraphWin(figName);
[y,x] = hist(nullDistribution,100);
y = y./sum(y);
bar(x,y,'k')
hold on
plot([EmpiricalDiff,EmpiricalDiff],[0 max(y)],'r-','linewidth',2)
set(gca,'tickdir','out','box','off', 'ylim',[0 max(y)],'FontSize',16)
ylabel('Probability')
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
saveFig(fh,fullfile(saveDir,figName),1)
 
% Make a plot of the R-squared
figName = sprintf('Test_ILF_SUBTRACTION_rmse_%s_lmax%i',fname,lmax(i_lmax));
fh = mrvNewGraphWin('Test ILF removal');
bar([WITH.rmsem,WITHOUT.rmsem],'FaceColor',colors{1})
hold on
plot([1 1; 2 2]',[WITH.rmsesd,WITHOUT.rmsesd],'r-','linewidth',16)
ylabel('rmse')
set(gca,'xticklabel',{'With fascicle','Without fascicle'},'tickdir','out','box','off', ...
    'ylim',[20 50],'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),1)

% Make a plot of the R-squared
figName = sprintf('Test_ILF_SUBTRACTION_Rrmse_%s_lmax%i',fname,lmax(i_lmax));
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
figName = sprintf('Test_ILF_SUBTRACTION_rmse_SCATTER_%s_lmax%i',fname,lmax(i_lmax));
fh = mrvNewGraphWin(figName);
hold on
plot([0 80],[0 80],'k-',[median(WITHOUT.rmse),median(WITHOUT.rmse)],[0 80],'k--',[0 80],[median(WITH.rmse),median(WITH.rmse)],'k--')

for ii = 1:2%length(WITHOUT.rmseall)
plot(WITHOUT.rmseall{ii},WITH.rmseall{ii},'o','color',colors{ii},'markerfacecolor',colors{ii});
end
axis square
set(gca,'tickdir','out','box','off','xlim',[0 80],'ylim',[0 80],'FontSize',16)
ylabel('rmse WITH ILF')
xlabel('rmse WITHOUT ILF')
saveFig(fh,fullfile(saveDir,figName),1)

end % lmax

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
