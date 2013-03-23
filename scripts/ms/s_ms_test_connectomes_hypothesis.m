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
if notDefined('trackingType'),trackingType = 'd';end
if notDefined('lmax'),        lmax         = 2;end
if notDefined('bval'),        bval         = 4000;end
if notDefined('rep'),         rep          = [1];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('cullType'),   cullType='culledL2';end

if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_figures');end
% High-resolution Anatomy
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';

% Path to the ROIs used to segment the fascicle
roisPath = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/ROIs_restric_connectomes';

% These are the ROIs used to segment the fascicle
rois = {fullfile(roisPath,'ROI_Anterior.mat'),  ...
    fullfile(roisPath,'ROI_Inferior.mat'),  ...
    fullfile(roisPath,'ROI_posterior.mat'),  ...
    fullfile(roisPath,'ROI_Superior.mat')};

% These are the logical operations to apply to each ROI
operation = {'and','not','and','not'};

% Make some plots
colors = {[.5 .4 .65],[.6 .45 .4],[.45 .6 .35]};

for irep = 1:length(rep)
    if notDefined('feFileToLoad')
        % Information on the path to the files to load.
        % This is where the inputs will be loaded from
        [feFileToLoad, fname] = msBuildFeFileName(trackingType,lmax,bval,rep(irep),diffusionModelParams,cullType);
    else
    fname = 'input';    
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
    
    % Perform a series of end and not operations to segement the OR/ILF
    [~, keepFascicles] = feSegmentFascicleFromConnectome(fg, rois, operation, 'prob connectome');
    
    % Remove all the voxels from the connectome except the ones where the
    % fascicle passes through. Fit the new model.
    [feWithoutFas, feWithFas, ~] = feTestFascicle(fe,keepFascicles,0);
    
    % Make a plot of the R-squared
    WITH.r2(irep)       = median(feGetRep(feWithFas,   'vox  r2'));
    WITH.rmse(irep)     = median(feGetRep(feWithFas,   'vox  rmse'));
    WITH.rrmse(irep)    = median(feGetRep(feWithFas,   'vox  rmse ratio'));
    
    WITHOUT.r2(irep)    = median(feGetRep(feWithoutFas,'vox  r2'));
    WITHOUT.rmse(irep)  = median(feGetRep(feWithoutFas,'vox  rmse'));
    WITHOUT.rrmse(irep) = median(feGetRep(feWithoutFas,'vox  rmse ratio'));
    
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

% Make figures and save them to file
figName = sprintf('Test_ILF_removal_r2_%s',fname);
fh = mrvNewGraphWin(figName);
bar(100.*[WITH.r2m,WITHOUT.r2m],'FaceColor',colors{1})
hold on
plot([1 1; 2 2]',100.*[WITH.r2sd,WITHOUT.r2sd],'r-','linewidth',16)
ylabel('Percent variance expalined')
set(gca,'xticklabel',{'With fascicle','Without fascicle'},'tickdir','out','box','off', ...
    'ylim',[35 60],'FontSize',16)
saveFig(fh,fullfile(saveDir,figName))
 
% Make a plot of the R-squared
figName = sprintf('Test_ILF_removal_rmse_%s',fname);
fh = mrvNewGraphWin('Test ILF removal');
bar([WITH.rmsem,WITHOUT.rmsem],'FaceColor',colors{1})
hold on
plot([1 1; 2 2]',[WITH.rmsesd,WITHOUT.rmsesd],'r-','linewidth',16)
ylabel('rmse')
set(gca,'xticklabel',{'With fascicle','Without fascicle'},'tickdir','out','box','off', ...
    'ylim',[20 40],'FontSize',16)
saveFig(fh,fullfile(saveDir,figName))

% Make a plot of the R-squared
figName = sprintf('Test_ILF_removal_Rrmse_%s',fname);
fh = mrvNewGraphWin('Test ILF removal');
bar([WITH.rrmsem,WITHOUT.rrmsem],'FaceColor',colors{1})
hold on
plot([1 1; 2 2]',[WITH.rrmsesd,WITHOUT.rrmsesd],'r-','linewidth',16)
plot([0 3],[1 1],'k-')
ylabel('R_{rmse}')
set(gca,'xticklabel',{'With fascicle','Without fascicle'},'tickdir','out','box','off', ...
    'ylim',[0.75 1.25],'FontSize',16)
saveFig(fh,fullfile(saveDir,figName))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
eval( sprintf('print(%s, ''-depsc2'',''-tiff'', ''%s'');', num2str(h),figName));

end
