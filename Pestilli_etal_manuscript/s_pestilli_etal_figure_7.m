function s_pestilli_etal_figure_7()
%% Comparison between two tractography model.
%
%  This function illustrates how to:
%  - Compute the Room-Mean-Square-Error (RMSE) from precomputed LiFE structures.
%  - Compare the RMSE of two different tractography models. In this example
%    we show how to compare a Probabilistic and a Deterministic
%    tractogrpahy model.
%  - We show how to use the LiFE software to compute the "Strength of evidence" and
%    the "Earth Movers Distance" to assess the difference in RMSE between
%    two tractography models.
%
% The example shows that for this brain and tractography model the
% Probabilitic tractorgrpahy model generates smaller error; it predicts the
% diffusion measuments better than the Deterministic.
% 
% Results similar to the ones in reproduced in this exmple are reported in
% Fig. 7 of Pestilli et al. UNDER REVIEW.
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com

% Get the base directory for the data
datapath = pestilliDataPath;

%% Load two pre-computed connectomes generated using Probabilistic and Deterministic tractography.
%
feProbFileName = 'subject1_life_culled_2mm_150dir_b2000_probabilistic_lmax8_diffModAx100Rd0.mat';
feDetFileName  = 'subject1_life_culled_2mm_150dir_b2000_tensor_diffModAx100Rd0.mat';

fprintf('Loading precomputed LiFE models for probabilistic (P) and deterministic (D) connectomes...\n')
p = load(fullfile(datapath,'life_structures',feProbFileName));
d = load(fullfile(datapath,'life_structures',feDetFileName)); 

%% Extract the RMSE of the LiFE model for each connectome in each white-matter voxel.
% We compute the RMSE for each individual voxel within the White-matter.
fprintf('Computing Root-Mean-Square-Error (RMSE) for each brain voxel and tractography model..\n')
p.rmse   = feGetRep(p.fe,'vox rmse');
d.rmse   = feGetRep(d.fe, 'vox rmse');

%% Find the common coordinates between the two connectomes.
% The two tractography method might have passed through slightly different
% voxels. Here we find the voxels where both models passed. We will compare
% the error only in these common voxels. There are more coordinates in the
% Prob connectome, because the tracking fills up more White-matter. 
%
% So, hereafter:
%
% - First we find the indices in the probabilistic connectome of the
% coordinate in the deterministic connectome. But there are some of the
% coordinates in the Deterministic conectome that are NOT in the
% Probabilistic connectome.
%
% - Second we find the indices in the Deterministic connectome of the
% subset of coordinates in the Probabilistic connectome found in the
% previous step.
%
% - Third we find the common voxels. These allow us to find the rmse for
% the same voxels.

fprintf('Finding common brain coordinates between P and D connectomes...\n')
p.coords = feGet(   p.fe,'roi coords');
d.coords = feGet(   d.fe, 'roi coords');
prob.coordsIdx = ismember(p.coords,d.coords,'rows');
prob.coords   = p.coords(prob.coordsIdx,:);
det.coordsIdx = ismember(d.coords,prob.coords,'rows');
det.coords    = d.coords(det.coordsIdx,:);
prob.rmse  = p.rmse( prob.coordsIdx);
det.rmse   = d.rmse( det.coordsIdx);


%% Compare the RMSE of the Probabilistic and Deterministic models using a scatter-density plot. 
scatterPlotRMSE(det,prob)

%% Compare the RMSE of the two models using the Stregth-of-evidence and the Earth Movers Distance.
se = feComputeEvidence(p.rmse,d.rmse);

%% Show the strength of evidence in favor of Probabilistic versus Deterministic tractography. 
% Plot the distributions of resampled mean RMSE used to compute the strength of
% evidence (S). 
distributionPlotStrengthOfEvidence(se.s.unlesioned_e,se.s.lesioned_e,se.s.mean,se.s.std,se.s.unlesioned.xbins, se.s.lesioned.xbins)

%% Show the RMSE distributions for Probabilistic deterministic tractography. 
%% Compare the distributions using the Earth Movers Distance.  
% Plot the distributions of RMSE for the two models and report the Earth
% Movers Distance between the distributions

distributionPlotEarthMoversDistance(se.nolesion,se.lesion,se.em)

end % End main function

%% Helper functions used for plotting
function scatterPlotRMSE(det,prob)
figNameRmse = sprintf('prob_vs_det_rmse_common_voxels_map');
mrvNewGraphWin(figNameRmse);
[ymap,x]  = hist3([det.rmse;prob.rmse]',{[10:1:70], [10:1:70]});
ymap = ymap./length(prob.rmse);
sh   = imagesc(flipud(log10(ymap)));
cm   = colormap(flipud(hot)); view(0,90);
axis('square')      
set(gca, ...
    'xlim',[1 length(x{1})],...
    'ylim',[1 length(x{1})], ...
    'ytick',[1 (length(x{1})/2) length(x{1})], ...
    'xtick',[1 (length(x{1})/2) length(x{1})], ...
    'yticklabel',[x{1}(end) x{1}(round(end/2)) x{1}(1)], ...
    'xticklabel',[x{1}(1)   x{1}(round(end/2)) x{1}(end)], ...
    'tickdir','out','ticklen',[.025 .05],'box','off', ...
    'fontsize',16,'visible','on')
hold on
plot3([1 length(x{1})],[length(x{1}) 1],[max(ymap(:)) max(ymap(:))],'k-','linewidth',1)
ylabel('Deterministic_{rmse}','fontsize',16)
xlabel('Probabilistic_{rmse}','fontsize',16)
cb = colorbar;
tck = get(cb,'ytick');
set(cb,'yTick',[min(tck)  mean(tck) max(tck)], ...
    'yTickLabel',round(1000*10.^[min(tck),...
    mean(tck), ...
    max(tck)])/1000, ...
    'tickdir','out','ticklen',[.025 .05],'box','on', ...
    'fontsize',16,'visible','on')
end

function distributionPlotStrengthOfEvidence(y_e,ywo_e,dprime,std_dprime,xhis,woxhis)
histcolor{1} = [0 0 0];
histcolor{2} = [.95 .6 .5];
figName = sprintf('Strength_of_Evidence_test_PROB_vs_DET_model_rmse_mean_HIST');
mrvNewGraphWin(figName);
patch([xhis,xhis],y_e(:),histcolor{1},'FaceColor',histcolor{1},'EdgeColor',histcolor{1});
hold on
patch([woxhis,woxhis],ywo_e(:),histcolor{2},'FaceColor',histcolor{2},'EdgeColor',histcolor{2}); 
set(gca,'tickdir','out', ...
        'box','off', ...
        'ticklen',[.025 .05], ...
        'ylim',[0 .2], ... 
        'xlim',[28 34], ...
        'xtick',[28 30 32 34], ...
        'ytick',[0 .1 .2], ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('rmse','fontsize',16')

title(sprintf('Strength of evidence:\n mean %2.3f - std %2.3f',dprime,std_dprime), ...
    'FontSize',16)
legend({'Probabilistic','Deterministic'})
end


function distributionPlotEarthMoversDistance(prob,det,em)
histcolor{1} = [0 0 0];
histcolor{2} = [.95 .6 .5];
figName = sprintf('EMD_PROB_DET_model_rmse_mean_HIST');
mrvNewGraphWin(figName);
plot(prob.xhist,prob.hist,'r-','color',histcolor{1},'linewidth',4);
hold on
plot(det.xhist,det.hist,'r-','color',histcolor{2},'linewidth',4); 
set(gca,'tickdir','out', ...
        'box','off', ...
        'ticklen',[.025 .05], ...
        'ylim',[0 .12], ... 
        'xlim',[0 95], ...
        'xtick',[0 45 90], ...
        'ytick',[0 .06 .12], ...
        'fontsize',16)
ylabel('Proportion white-matter volume','fontsize',16)
xlabel('RMSE (raw MRI scanner units)','fontsize',16')
title(sprintf('Earth Movers Distance: %2.3f (raw scanner units)',em.mean),'FontSize',16)
legend({'Probabilistic','Deterministic'})
end