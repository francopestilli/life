function s_ms_test_fiber_angles_TEMP_procomputed(feFiles,saveDir)
%
% s_ms_test_fiber_angles(feFiles,saveDir)
%
% Load a series of FE structures for conectomes of the occipital lobe.
% Plots the distribution of angles across fiebrs in the conenctome.
%
% This is part of a series of reproducible science scritps to be published
% with the LiFE mansucript.
%
% Franco (c) Stanford Vista Team 2013
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_figures');end

% Make some plots
colors = {[.9 .45 .25], [.5 .4 .65],[.45 .6 .35]};

% Initalize the number the vectors that will contain the distribution fo
% angles for the connectome.
nBars = 180; % Number of bins in the histogram plots.
ely = nan(1,nBars);
elx = nan(nBars,1);
azy = ely;azx = elx;
for irep = 1:length(feFiles)
    % Extract the fiber group from the FE structure
    fg = feGet(feFiles{irep},'fibers img');
    [~,~,azv,elv,azvx,elvx] = feComputeFiberAngle(fg);
    
    % Now compute the distribution fo angles across trakings.
    [ely(irep,:),elx] = hist(elv,nBars);
    [azy(irep,:),azx] = hist(azv,nBars);
end

% Compute the entropy for each histogram:
% To do so we transfrom the hisotgram in probabilities.
% Elevation.
tot_el     =  sum(ely,2);
ely_p      =  ely./repmat(tot_el,1,size(ely,2));
entropy_el = -sum(ely_p.*log2(ely_p),2);

% Azimuth.
tot_az     =  sum(azy,2);
azy_p      =  azy./repmat(tot_az,1,size(azy,2));
entropy_az = -sum(azy_p.*log2(azy_p),2);
    
% Make figures and save them to file
for irep = 1:length(feFiles)
    if irep == 1
        figNameE = sprintf('Elevation_Prob_Det_Tensor');
        fhE = mrvNewGraphWin(figNameE);
    else
        figure(fhE)
    end
    hold on
    hbe(irep) = bar(elx,ely(irep,:)./max(ely(irep,:)), ...
        'FaceColor',colors{irep},'EdgeColor',colors{irep});
    set(get(hbe(irep),'Children'),'FaceAlpha',.5,'EdgeAlpha',.5);
    ylabel('Frequency of occurrence (normalized)');
    xlabel('Elevation (degrees)');
    set(gca,'tickdir','out','box','off', 'FontSize',16);
end
legend(hbe,{sprintf('Probabilistic (entropy %2.3f)',entropy_el(1)), ...
            sprintf('Deterministic (entropy %2.3f)',entropy_el(2)), ...
            sprintf('Tensor        (entropy %2.3f)',entropy_el(3))},'box','off')
saveFig(fhE,fullfile(saveDir,figNameE));
    
% Make a plot of the R-squared
for irep = 1:length(feFiles)
    if irep == 1
        figNameA = sprintf('Azimuth_Prob_Det_Tensor');
        fhA = mrvNewGraphWin(figNameA);
    else
        figure(fhA)
    end
    hold on
    hba(irep) = bar(azx,azy(irep,:)./max(azy(irep,:)), ...
         'FaceColor', colors{irep},'EdgeColor',colors{irep});
    set(get(hba(irep),'Children'),'FaceAlpha',.5,'EdgeAlpha',.5);
    ylabel('Frequency of occurrence (normalized)');
    xlabel('Azimuth (degrees)');
    set(gca,'tickdir','out','box','off','FontSize',16);
end
legend(hba,{sprintf('Probabilistic (entropy %2.3f)',entropy_az(1)), ...
            sprintf('Deterministic (entropy %2.3f)',entropy_az(2)), ...
            sprintf('Tensor        (entropy %2.3f)',entropy_az(3))},'box','off')
saveFig(fhA,fullfile(saveDir,figNameA));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
eval( sprintf('print(%s,  ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'');', num2str(h),figName));

end
