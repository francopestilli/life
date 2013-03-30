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
if notDefined('trackingType'),trackingType = 't';end
if notDefined('lmax'),        lmax         = 6;end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = [1,2,3];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('cullType'),   cullType='culledL2';end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_fiber_angles');end
if notDefined('feFileToLoad'), feFileToLoad = 'load from disk';
else rep =1;
end
% Check that the connectome were proprocessed before attempting to make a
% plot.
done = s_ms_check_processes([],trackingType,lmax,bval,cullType);
if ~all(done), 
    fprintf('\n[%s] Not all connectomes were proprocessed... not proceeding with plot.\n',mfilename);
return
end

% Make some plots
colors = {[.9 .45 .25], [.5 .4 .65],[.45 .6 .35]};

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
    [~,~,azv,elv] = feComputeFiberAngle(fg);
    
    % Now compute the distribution fo angles across trakings.
    [ely(irep,:),elx] = hist(elv,nBars);
    [azy(irep,:),azx] = hist(azv,nBars);
end

% Elevation.
tot_el     =  sum(ely,2);
ely_p      =  ely./repmat(tot_el,1,size(ely,2));
entropy_el = -sum(ely_p.*log2(ely_p),2);

% Azimuth.
tot_az     =  sum(azy,2);
azy_p      =  azy./repmat(tot_az,1,size(azy,2));
entropy_az = -sum(azy_p.*log2(azy_p),2);

% Make entropy plots:
figName = sprintf('Elevation_entropy_%s',fname);
fh = mrvNewGraphWin(figName);
bar(mean(entropy_el),'FaceColor',colors{1},'EdgeColor',colors{1});
hold on
plot([1;1],[mean(entropy_el); ...
     mean(entropy_el)] +[-std(entropy_el); std(entropy_el)], ...
     'k-','linewidth',10)
set(gca,'xlim',[0,2],'ylim',[6.75,7.75],'tickdir','out','box','off', 'FontSize',16);
saveFig(fh,fullfile(saveDir,figName));

figName = sprintf('Azimuth_entropy_%s',fname);
fh = mrvNewGraphWin(figName);
bar(mean(entropy_az),'FaceColor',colors{1},'EdgeColor',colors{1});
hold on
plot([1;1],[mean(entropy_az); ...
     mean(entropy_az)] +[-std(entropy_az); std(entropy_az)], ...
     'k-','linewidth',10)
set(gca,'xlim',[0,2],'ylim',[6.75,7.75],'tickdir','out','box','off', 'FontSize',16);
saveFig(fh,fullfile(saveDir,figName));

% Make figures and save them to file
figName = sprintf('Elevation_%s',fname);
fh = mrvNewGraphWin(figName);
hb = bar(elx,mean(ely./size(azy,2),1),'FaceColor',colors{1},'EdgeColor',colors{1});
set(get(hb,'Children'),'FaceAlpha',.5,'EdgeAlpha',.5);
hold on
plot([elx;elx],[mean(ely./size(azy,2),1); ...
     mean(ely./size(azy,2),1)] +[-std(ely./size(azy,2),[],1); std(ely./size(azy,2),[],1)], ...
     'r-','linewidth',2)
ylabel('Probability of angle');
xlabel('Elevation in Degrees');
set(gca,'tickdir','out','box','off', ...
    'FontSize',16);
legend(hb,{sprintf('Entropy %2.3f (sd = %2.3f)',mean(entropy_el),std(entropy_el))},'box','off')
saveFig(fh,fullfile(saveDir,figName));
 
% Make a plot of the R-squared
figName = sprintf('Azimuth_%s',fname);
fh = mrvNewGraphWin(figName);
hb = bar(azx,mean(azy./size(azy,2),1),'FaceColor',colors{1},'EdgeColor',colors{1});
set(get(hb,'Children'),'FaceAlpha',.5,'EdgeAlpha',.5);
hold on
plot([azx;azx],[mean(azy./size(azy,2),1); ...
     mean(azy./size(azy,2),1)] +[-std(azy./size(azy,2),[],1); std(azy./size(azy,2),[],1)], ...
     'r-','linewidth',2)
ylabel('Probability of angle');
xlabel('Azimuth in Degrees');
set(gca,'tickdir','out','box','off', ...
    'FontSize',16);
legend(hb,{sprintf('Entropy %2.3f (sd = %2.3f)',mean(entropy_az),std(entropy_az))},'box','off')
saveFig(fh,fullfile(saveDir,figName));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end
eval( sprintf('print(%s,  ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'');', num2str(h),figName));

end
