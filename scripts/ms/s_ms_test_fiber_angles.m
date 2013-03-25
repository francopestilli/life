function s_ms_test_fiber_angles(feFileToLoad,trackingType,lmax,bval,rep)
%
% s_ms_test_fiber_angles(trackingType,lmax,bval,rep)
%
% Load a series of FE structures for conectomes of the occipital lobe.
% Plots the distribution of angles across fiebrs in the conenctome.
%
% This is part of a series of reproducible science scritps to be published
% with the LiFE mansucript.
%
% Franco (c) Stanford Vista Team 2013
if notDefined('trackingType'),trackingType = 'p';end
if notDefined('lmax'),        lmax         = 2;end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = [1,2,3];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('cullType'),   cullType='culledL2';end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_figures');end
if notDefined('feFileToLoad'), feFileToLoad = 'load from disk';end

% Make some plots
colors = {[.6 .45 .4],[.5 .4 .65],[.45 .6 .35]};

% Initalize the number the vectors that will contain the distribution fo
% angles for the connectome.
nBars = 180; % Number of bins in the histogram plots.
ely = nan(length(rep),nBars);
elx = nan(length(rep),1);
azy = ely;azx = elx;
for irep = 1:length(rep)
    % Get the fe structure
    if ischar(feFileToLoad)
            % Information on the path to the files to load.
    % This is where the inputs will be loaded from
    [feFileToLoad, fname] = ...
        msBuildFeFileName(trackingType,lmax,bval,rep(irep), ...
        diffusionModelParams,cullType);
    fprintf('Loading %s ...\n',feFileToLoad)
        
        load(feFileToLoad);
    else
        fe  =feFileToLoad;
        clear feFileToLoad;
        fname = fe.name;
    end
    
    % Extract the fiber group from the FE structure
    fg = feGet(fe,'fibers acpc');
    [~,~,azv,elv] = feComputeFiberAngle(fg.fibers);
    
    % Now compute the distribution fo angles across trakings.
    [ely(irep,:),elx] = hist(elv,nBars);
    [azy(irep,:),azx] = hist(azv,nBars);
end

% Make figures and save them to file
figName = sprintf('Elevation_%s',fname);
fh = mrvNewGraphWin(figName);
bar(elx,mean(ely./size(azy,2)),'FaceColor',colors{1},'EdgeColor','w');
hold on
plot([elx;elx],[mean(ely./size(azy,2)); ...
     mean(ely./size(azy,2))] +[-std(ely./size(azy,2)); std(ely./size(azy,2))], ...
     'r-','linewidth',3)
ylabel('Probability of angle');
xlabel('Elevation in Degrees');
set(gca,'tickdir','out','box','off', ...
    'FontSize',16);
saveFig(fh,fullfile(saveDir,figName));
 
% Make a plot of the R-squared
figName = sprintf('Azimuth_%s',fname);
fh = mrvNewGraphWin(figName);
bar(azx,mean(azy./size(azy,2)),'FaceColor',colors{1},'EdgeColor','w');
hold on
plot([azx;azx],[mean(azy./size(azy,2)); ...
     mean(azy./size(azy,2))] +[-std(azy./size(azy,2)); std(azy./size(azy,2))], ...
     'r-','linewidth',3)
ylabel('Probability of angle');
xlabel('Azimuth in Degrees');
set(gca,'tickdir','out','box','off', ...
    'FontSize',16);
saveFig(fh,fullfile(saveDir,figName));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
eval( sprintf('print(%s, ''-depsc2'',''-tiff'', ''%s'');', num2str(h),figName));

end
