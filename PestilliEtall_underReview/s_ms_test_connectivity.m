function s_ms_test_connectivity(feFileToLoad,trackingType,lmax,bval,rep)
%
% s_ms_test_connectivity(trackingType,lmax,bval,rep)
%
% Load a series of FE structures for conectomes of the occipital lobe.
% Selectes the fibers touching all V1 and V2 voxesl.
%
% Computes a measure of conenctivity:
%
% the sum of the weights of the fibers touching the ROI, divided by the sum
% of the weights of all the other fibers going throught all the voxels of
% the fibers touching the ROI.
%
% This is part of a series of reproducible science scritps to be published
% with the LiFE mansucript.
%
% Franco (c) Stanford Vista Team 2013
if notDefined('trackingType'),trackingType = 't';end
if notDefined('lmax'),        lmax         = 2;end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = [1,2,3];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('cullType'),   cullType='culledL2';end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_figures');end

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
    % Information on the path to the files to load.
    % This is where the inputs will be loaded from
    [feFileToLoad, fname] = msBuildFeFileName(trackingType,lmax,bval,rep(irep),diffusionModelParams,cullType);
    
    % Get the fe structure
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
    if all(keepFascicles==0)
        keyboard
    end
    
    % Remove all the voxels from the connectome except the ones where the
    % fascicle passes through. Fit the new model.
    [~, ~, connectivity(irep)] = feTestFascicle(fe,keepFascicles,0);
end

% Now make averages and std or the connectivity estimates
for irep = 1:length(rep)
    wRatio(irep) = connectivity(irep).wratio;
    varExp(irep) = connectivity(irep).varexp;
end

% Compute mean and standard deviation
meanRw = mean(wRatio);
errRw  = [meanRw - std(wRatio), meanRw + std(wRatio)];
meanVE = mean(varExp);
errVE  = [meanVE - std(varExp), meanVE + std(varExp)];

% Make figures and save them to file
figName = sprintf('Connectivity_VarianceExplained_%s',fname);
fh = mrvNewGraphWin(figName);
bar(meanVE,'FaceColor',colors{1})
hold on
plot([1 1],errVE,'r-','linewidth',16)
ylabel('Percent variance expalined')
set(gca,'xticklabel',{'With fascicle','Without fascicle'},'tickdir','out','box','off', ...
    'ylim',[0 60],'FontSize',16)
saveFig(fh,fullfile(saveDir,figName))
 
% Make a plot of the R-squared
figName = sprintf('Connectivity_FascicleContributionRatio_%s',fname);
fh = mrvNewGraphWin(figName);
bar(meanRw,'FaceColor',colors{1})
hold on
plot([1 1],errRw,'r-','linewidth',16)
ylabel('R_{fascicle contribution}')
set(gca,'xticklabel',{'With fascicle','Without fascicle'},'tickdir','out','box','off', ...
    'ylim',[0 4],'ytick',[0 1 2 4],'FontSize',16,'yscale','lin')
saveFig(fh,fullfile(saveDir,figName))


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
eval( sprintf('print(%s, ''-depsc2'',''-tiff'', ''%s'');', num2str(h),figName));

end
