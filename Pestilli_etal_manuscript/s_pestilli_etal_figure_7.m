function s_pestilli_etal_figure_7()
%
% This function llustrates how to:
%  - initialize a LIFE structure from a candidate connectome
%  - Generate an optimized connectome from a cadidate connectome using the
%  LIFE strustrue
%
%  fe = s_pestilli_etal_figure_7()
% 
% INPUTS:  none
% OUTPUTS: fe structure the optimized life structure
%
% Copyright 2013-2014 Franco Pestilli Stanford University pestillifranco@gmail.com.

% Get the base directory for the data
datapath = pestilliDataPath;


%% Load two pre-culled connectomes
%
feProbFileName = 'subject1_life_culled_2mm_150dir_b2000_probabilistic_lmax8_diffModAx100Rd0.mat';
feDetFileName  = 'subject1_life_culled_2mm_150dir_b2000_tensor_diffModAx100Rd0.mat';

p = load(fullfile(datapath,'life_structures',feProbFileName));
d = load(fullfile(datapath,'life_structures',feDetFileName)); 

% Extract the RMSE, R_rmse and the coordinates of the white matter
p.rmse   = feGetRep(p.fe,'vox rmse');
p.rrmse  = feGetRep(p.fe,'vox rmse ratio');
p.coords = feGet(   p.fe,'roi coords');

d.rmse   = feGetRep(d.fe, 'vox rmse');
d.rrmse  = feGetRep(d.fe, 'vox rmse ratio');
d.coords = feGet(   d.fe, 'roi coords');


%% Find the common coordinates between the two connectomes
%
% There are more coordinates in the Prob conectome, because the tracking
% fills up more White-matter.
%
% So, first we find the indices in the probabilistic connectome of the
% coordinate in the deterministic conenctome.
%
% But there are some of the coordinates in the Deterministic conectome that
% are NOT in the Probabilistic connectome.
% 
% So, second we find the indices in the Deterministic connectome of the
% subset of coordinates in the Probabilistic connectome found in the
% previous step.

% First we find the coordinates in the Probabilistic conectome that are
% also in the Deterministic connectome.
prob.coordsIdx = ismember(p.coords,d.coords,'rows');

% Second we find the coordinates in the Deterministic connectome that are
% also in the Probabilistic connectome.
prob.coords   = p.coords(prob.coordsIdx,:);
det.coordsIdx = ismember(d.coords,prob.coords,'rows');
det.coords    = d.coords(det.coordsIdx,:);

% What we really need is detCoordsIdx and probCoordsIdx. These allow us to
% find the common voxel indices in rmse and rrmse, etc.
prob.rmse  = p.rmse( prob.coordsIdx);
prob.rrmse = p.rrmse(prob.coordsIdx);

det.rmse  = d.rmse( det.coordsIdx);
det.rrmse = d.rrmse(det.coordsIdx);


%% RMSE scatter-density plot of  Probabilistic and Deterministic
scatterPlotRMSE(det,prob)


%% Make a statisitcal test. To show that the Probabilistic model is better
% than the deterministic model when using all the voxels.
nmontecarlo = 5;
nboots      = 1000;
nbins       = 200;
sizeWith    = length(det.rmse);
nullDistributionP = nan(nboots,nmontecarlo);
nullDistributionD = nan(nboots,nmontecarlo);
y = nan(nbins,nmontecarlo);woy = y;

% Repeat the bootstrap several times
for inm = 1:nmontecarlo
    % Bootstrap the mean RMSE
    parfor ibt = 1:nboots
        nullDistributionP(ibt,inm) = mean(randsample(prob.rmse, sizeWith,true));      
        nullDistributionD(ibt,inm) = mean(randsample(det.rmse, sizeWith,true));
    end
    
    % Distribution probabilistic
    [y(:,inm),xhis] = hist(nullDistributionP(:,inm),linspace(22,34,200));
    y(:,inm) = y(:,inm)./sum(y(:,inm));
    
    % Distribution deterministic
    [woy(:,inm),woxhis] = hist(nullDistributionD(:,inm),linspace(22,34,200));
    woy(:,inm) = woy(:,inm)./sum(woy(:,inm));
end
y_m = mean(y,2);
y_e = [y_m, y_m] + 2*[-std(y,[],2),std(y,[],2)];

ywo_m = mean(woy,2);
ywo_e = [ywo_m, ywo_m] + 2*[-std(woy,[],2),std(woy,[],2)];


%% Compute the strength of evidence.
dprime = diff([mean(nullDistributionP,1);mean(nullDistributionD,1)]) ...
               ./sqrt(sum([std(nullDistributionP,[],1);std(nullDistributionD,[],1)].^2,1));

           
%% Plot the null distribution and the empirical difference
distributionPlotRMSE(y_e,ywo_e,dprime,xhis,woxhis)


end

%------------------------------------%
function scatterPlotRMSE(det,prob)
figNameRmse = sprintf('prob_vs_det_rmse_common_voxels_map');
fhRmseMap = mrvNewGraphWin(figNameRmse);
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
ylabel('Deterministic_{rmse}','fontsize',12)
xlabel('Probabilistic_{rmse}','fontsize',12)
cb = colorbar;
tck = get(cb,'ytick');
set(cb,'yTick',[min(tck)  mean(tck) max(tck)], ...
    'yTickLabel',round(1000*10.^[min(tck),...
    mean(tck), ...
    max(tck)])/1000, ...
    'tickdir','out','ticklen',[.025 .05],'box','on', ...
    'fontsize',16,'visible','on')
end

%---------------------------------------%
function distributionPlotRMSE(y_e,ywo_e,dprime,xhis,woxhis)

h1.ylim  = [0 0.6];
h1.xlim  = [22,34];
h1.ytick = [0 0.3 0.6];
h1.xtick = [28 30 32 34];
h2.ylim  = [0 0.4];
h2.xlim  = [28,32];
h2.ytick = [0 0.2 0.4];
h2.xtick = [28 30 32];
histcolor{1} = [0 0 0];
histcolor{2} = [.95 .6 .5];

figName = sprintf('Test_PROB_DET_model_rmse_mean_HIST');
fh = mrvNewGraphWin(figName);
patch([xhis,xhis],y_e(:),histcolor{1},'FaceColor',histcolor{1},'EdgeColor',histcolor{1});
hold on
patch([woxhis,woxhis],ywo_e(:),histcolor{2},'FaceColor',histcolor{2},'EdgeColor',histcolor{2}); 
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[0 .6], ... 
        'xlim',[28 34], ...
        'xtick',[28 30 32 34], ...
        'ytick',[0 .3 .6], ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('rmse','fontsize',16')

title(sprintf('Strength of evidence:\n mean %2.3f - std %2.3f',mean(dprime),std(dprime)), ...
    'FontSize',16)
legend({'Probabilistic','Deterministic'})
end