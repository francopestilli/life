function feFileToLoad = s_ms_connectomes_rmse_compare_tractography(feType,connectomeType,rep)
%
% feFileToSave = s_ms_connectomes_rmse_compare_tractography(connectomeType,rep)
%
% Loads a connectome and produces a figure of the rmse M-D and rmse D-D.
%
% Written by Franco Pestilli (c) Stanford University 2013 
if notDefined('connectomeType'), connectomeType = [41 4];end
if notDefined('rep'),         rep          = [1];end
if notDefined('cullType'),   cullType='culledL2';end
if notDefined('feType'),feType = 'rh';end
hemisphere_path = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/fe_hemispheres';
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_rmse_compare',feType);end

fontSiz = 16;
nboots = 10000;
nmontecarlo = 10;
histcolor{1} = [0 0 0];
histcolor{2} = [.95 .6 .5];

for ic = 1:length(connectomeType)
    switch feType
        case {'occipital','o'}
            [trackingType, lmax, bval, diffusionModelParams] = getConditions(connectomeType(ic));
            [feFileToLoad, fName] = msBuildFeFileName(trackingType,lmax,bval,rep,diffusionModelParams,cullType);
            bymax = 14;
            h1.ylim  = [0 0.25];
            h1.xlim  = [29,33];
            h1.ytick = [0 0.1 0.2];
            h1.xtick = [29 31 33];
            h2.ylim  = [0 0.25];
            h2.xlim  = [26,29];
            h2.ytick = [0 0.1 0.2];
            h2.xtick = [26 27 28 29];
            
        case {'lhprobtensor'}
            tensorProb{1} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_TENSOR_left.mat');
            tensorProb{2} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_LMAX8_left.mat');
            feFileToLoad = tensorProb{ic};
            bymax = 100;  
            h1.ylim  = [0 0.6];
            h1.xlim  = [28,34];
            h1.ytick = [0 0.3 0.6];
            h1.xtick = [28 30 32 34];
            h2.ylim  = [0 0.4];
            h2.xlim  = [28,32];
            h2.ytick = [0 0.2 0.4];
            h2.xtick = [28 30 32];
            
        case {'rhprobtensor'}
            tensorProb{1} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_TENSOR_right.mat');
            tensorProb{2} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_LMAX8_right.mat');
            feFileToLoad = tensorProb{ic};
            bymax = 100;
            h1.ylim  = [0 0.6];
            h1.xlim  = [28,34];
            h1.ytick = [0 0.3 0.6];
            h1.xtick = [28 30 32 34];
            h2.ylim  = [0 0.4];
            h2.xlim  = [28,32];
            h2.ytick = [0 0.2 0.4];
            h2.xtick = [28 30 32];
            
        case {'lhprobdet'}
            tensorProb{1} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_DET_LMAX2_left.mat');
            tensorProb{2} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_LMAX8_left.mat');
            feFileToLoad = tensorProb{ic};
            bymax = 100;
            h1.ylim  = [0 0.6];
            h1.xlim  = [26,36];
            h1.ytick = [0 0.3 0.6];
            h1.xtick = [26 31 36];
            h2.ylim  = [0 0.4];
            h2.xlim  = [28,32];
            h2.ytick = [0 0.2 0.4];
            h2.xtick = [28 30 32];
            
        case {'rhprobdet'}
            tensorProb{1} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_DET_LMAX2_right.mat');
            tensorProb{2} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_LMAX8_right.mat');
            feFileToLoad = tensorProb{ic};
            bymax = 100;  
            h1.ylim  = [0 0.6];
            h1.xlim  = [26,36];
            h1.ytick = [0 0.3 0.6];
            h1.xtick = [26 31 36];
            h2.ylim  = [0 0.4];
            h2.xlim  = [28,32];
            h2.ytick = [0 0.2 0.4];
            h2.xtick = [28 30 32];
            
        case {'lhtensdet'}
            tensorProb{1} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_TENSOR_left.mat');
            tensorProb{2} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_DET_LMAX2_left.mat');
            feFileToLoad = tensorProb{ic};
            bymax = 100;
            h1.ylim  = [0 0.6];
            h1.xlim  = [26,36];
            h1.ytick = [0 0.3 0.6];
            h1.xtick = [26 31 36];
            h2.ylim  = [0 0.4];
            h2.xlim  = [28,32];
            h2.ytick = [0 0.2 0.4];
            h2.xtick = [28 30 32];
            
        case {'rhtensdet'}
            tensorProb{1} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_TENSOR_right.mat');
            tensorProb{2} = fullfile(hemisphere_path,'fe_culled_FP_150_B2000_DET_LMAX2_right.mat');
            feFileToLoad = tensorProb{ic};
            bymax = 100;  
            h1.ylim  = [0 0.6];
            h1.xlim  = [26,36];
            h1.ytick = [0 0.3 0.6];
            h1.xtick = [26 31 36];
            h2.ylim  = [0 0.4];
            h2.xlim  = [28,32];
            h2.ytick = [0 0.2 0.4];
            h2.xtick = [28 30 32];
        otherwise
            keyboard
    end
    
    if exist(feFileToLoad,'file')
        fprintf('[%s] Found connectome, loading it:\n%s\n',mfilename,feFileToLoad);
        load(feFileToLoad)
    else
        error('[%s] Cannot find the FE structure...\n%s\n',mfilename,feFileToLoad);
    end
    
    % Extract the data to dat reliability:
    coords{ic} = feGet(fe,'roi coords');
    Drmse{ic}  = feGetRep(fe, 'vox rmse data');
    Mrmse{ic}  = feGetRep(fe, 'vox rmse');
    Rrmse{ic}  = feGetRep(fe, 'vox rmse ratio');
    Psig{ic}   = feGet(fe, 'pSig fiber');
    nbvas{ic}  = feGet(fe,'nbvals');
    nvoxs{ic}  = feGet(fe,'nvoxels');
    
end

   
% Find the common coordinates between the two connectomes
%
% There are more coordinates in the Prob conectome, because the tracking
% fills up more White-matter volume.
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

% First wwe find the coordinates in the Probabilistic conectome that are
% also in the Deterministic connectome.
prob.coordsIdx = ismember(coords{2},coords{1},'rows');

% Second we find the coordinates in the Deterministic connectome that are
% also in the Probabilistic connectome.
prob.coords   = coords{2}(prob.coordsIdx,:);
det.coordsIdx = ismember(coords{1},prob.coords,'rows');
det.coords    = coords{1}(det.coordsIdx,:);

% What we rally need is detCoordsIdx and probCoordsIdx. These allow us to
% find the common voxel indices in rmse and Rrmse, etc.
prob.Drmse = Drmse{2}(prob.coordsIdx);
prob.Mrmse = Mrmse{2}(prob.coordsIdx);
prob.Rrmse = Rrmse{2}(prob.coordsIdx);

det.Drmse = Drmse{1}(det.coordsIdx);
det.Mrmse = Mrmse{1}(det.coordsIdx);
det.Rrmse = Rrmse{1}(det.coordsIdx);

figNameRmse = sprintf('prob_vs_det_rmse_common_voxels_map');
fhRmseMap = mrvNewGraphWin(figNameRmse);
[ymap,x]  = hist3([det.Mrmse;prob.Mrmse]',{[10:1:70], [10:1:70]});
ymap = ymap./length(prob.Mrmse);
sh = imagesc(flipud(log10(ymap)));
cm = colormap(flipud(hot)); view(0,90);
axis('square')      
set(gca, ...
    'xlim',[1 length(x{1})],...
    'ylim',[1 length(x{1})], ...
    'ytick',[1 (length(x{1})/2) length(x{1})], ...
    'xtick',[1 (length(x{1})/2) length(x{1})], ...
    'yticklabel',[x{1}(end) x{1}(round(end/2)) x{1}(1)], ...
    'xticklabel',[x{1}(1)   x{1}(round(end/2)) x{1}(end)], ...
    'tickdir','out','ticklen',[.025 .05],'box','off', ...
    'fontsize',fontSiz','visible','on')
hold on
plot3([1 length(x{1})],[length(x{1}) 1],[max(ymap(:)) max(ymap(:))],'k-','linewidth',1)
ylabel('Deterministic_{rmse}','fontsize',fontSiz)
xlabel('Probabilistic_{rmse}','fontsize',fontSiz)
cb = colorbar;
tck = get(cb,'ytick');
set(cb,'yTick',[min(tck)  mean(tck) max(tck)], ...
    'yTickLabel',round(1000*10.^[min(tck),...
    mean(tck), ...
    max(tck)])/1000, ...
    'tickdir','out','ticklen',[.025 .05],'box','on', ...
    'fontsize',fontSiz','visible','on')
saveFig(fhRmseMap,fullfile(saveDir,figNameRmse),'eps')
    
figNameRrmse = sprintf('prob_vs_det_Rrmse_map_common_voxels');
fhRmseMap = mrvNewGraphWin(figNameRrmse);
[ymap,x]  = hist3([det.Rrmse;prob.Rrmse]',{[0.5:.025:2], [0.5:.025:2]});
ymap = ymap./length(prob.Rrmse);
sh = imagesc(flipud(log10(ymap)));
cm = colormap(flipud(hot)); view(0,90);
axis('square')      
set(gca, ...
    'xlim',[1 length(x{1})],...
    'ylim',[1 length(x{1})+.1], ...
    'ytick',[1 ceil(length(x{1})/2) length(x{1})], ...
    'xtick',[1 ceil(length(x{1})/2) length(x{1})], ...
    'yticklabel',[x{1}(end) 1 x{1}(1)], ...
    'xticklabel',[x{1}(1)   1 x{1}(end)], ...
    'tickdir','out','ticklen',[.025 .05],'box','off', ...
    'fontsize',fontSiz','visible','on')
hold on
plot3([1 length(x{1})],[length(x{1}) 1],[max(ymap(:)) max(ymap(:))],'k-','linewidth',1)
ylabel('Deterministic_{R}','fontsize',fontSiz)
xlabel('Probabilistic_{R}','fontsize',fontSiz)
cb = colorbar;
tck = get(cb,'ytick');
set(cb,'yTick',[min(tck)  mean(tck) max(tck)], ...
    'yTickLabel',round(1000*10.^[min(tck),...
    mean(tck), ...
    max(tck)])/1000, ...
    'tickdir','out','ticklen',[.025 .05],'box','on', ...
    'fontsize',fontSiz','visible','on')
saveFig(fhRmseMap,fullfile(saveDir,figNameRrmse),'eps')

% So far we have made plots by looking at the common voxels.
% But the voxels that are not in the Determinstic model (which has less voxels)
% can also be used for the comparisons. The prediction of the signal in
% these voxels returned by the Deterministic connectome is 0 at all diffusion 
% directions (only isotropic diffusion). Hereafter I build the prediction 
% for these voxels. 

% We reshape the pSignal vector into voxelsXmeasurements
for ic  = 1:length(Rrmse)
    psig{ic} = reshape(Psig{ic},nbvas{ic},nvoxs{ic});
end

% Voxels in the Probabilistic connectome that are
% also in the Deterministic connectome.
prob.coordsIdx ;

% Voxels in the Deterministic connectome that are
% also in the Probabilistic connectome.
det.coordsIdx ;

% Make a new vector for the psig of the Deterministic. THis vector will
% contain zeros at each voxel in Probabilistic not contained in
% Deterministic. It will contain det.psig in the rest of the voxles.
det_psig = zeros(size(psig{2}));
det_psig(:,prob.coordsIdx) = psig{1}(:,det.coordsIdx);

% Now extract the measured diffusion signal the RMSE. NOTE: 
% The measured signal is the one from the Probabilistic connectome 
% (from the second data set).
dSig = feGetRep(fe,'dsigdemeanedvox');

% Compute the RMSE
det.rmse_prob = sqrt( mean((dSig - det_psig).^2,1) );
prob.rmse_all = feGetRep(fe,'voxrmse');

% Make a scatter plot.
figNameRmse = sprintf('prob_vs_det_rmse_map_all_voxels');
fhRmseMap = mrvNewGraphWin(figNameRmse);
[ymap,x]  = hist3([det.rmse_prob;prob.rmse_all]',{[10:1:70], [10:1:70]});
ymap = ymap./length(prob.Mrmse);
sh = imagesc(flipud(log10(ymap)));
cm = colormap(flipud(hot)); view(0,90);
axis('square')      
set(gca, ...
    'xlim',[1 length(x{1})],...
    'ylim',[1 length(x{1})], ...
    'ytick',[1 (length(x{1})/2) length(x{1})], ...
    'xtick',[1 (length(x{1})/2) length(x{1})], ...
    'yticklabel',[x{1}(end) x{1}(round(end/2)) x{1}(1)], ...
    'xticklabel',[x{1}(1)   x{1}(round(end/2)) x{1}(end)], ...
    'tickdir','out','ticklen',[.025 .05],'box','off', ...
    'fontsize',fontSiz','visible','on')
hold on
plot3([1 length(x{1})],[length(x{1}) 1],[max(ymap(:)) max(ymap(:))],'k-','linewidth',1)
ylabel('Deterministic_{rmse}','fontsize',fontSiz)
xlabel('Probabilistic_{rmse}','fontsize',fontSiz)
cb = colorbar;
tck = get(cb,'ytick');
set(cb,'yTick',[min(tck)  mean(tck) max(tck)], ...
    'yTickLabel',round(1000*10.^[min(tck),...
    mean(tck), ...
    max(tck)])/1000, ...
    'tickdir','out','ticklen',[.025 .05],'box','on', ...
    'fontsize',fontSiz','visible','on')
saveFig(fhRmseMap,fullfile(saveDir,figNameRmse),'eps')
 
% Make a statisitcal test. To show that the Probabilistic model is better
% than the deterministic model when using all the voxels.
sizeWith    = length(det.rmse_prob);
nullDistributionP = nan(nboots,nmontecarlo);
nullDistributionD = nan(nboots,nmontecarlo);

for inm = 1:nmontecarlo
    parfor ibt = 1:nboots
        nullDistributionP(ibt,inm) = mean(randsample(prob.rmse_all, sizeWith,true));      
        nullDistributionD(ibt,inm) = mean(randsample(det.rmse_prob, sizeWith,true));
    end
    
    % Distribution With
    [y(:,inm),xhis] = hist(nullDistributionP(:,inm),linspace(h1.xlim(1),h1.xlim(2),200));
    y(:,inm) = y(:,inm)./sum(y(:,inm));
    
    % Distribution without
    [woy(:,inm),woxhis] = hist(nullDistributionD(:,inm),linspace(h1.xlim(1),h1.xlim(2),200));
    woy(:,inm) = woy(:,inm)./sum(woy(:,inm));
end
y_m = mean(y,2);
y_e = [y_m, y_m] + 2*[-std(y,[],2),std(y,[],2)];

ywo_m = mean(woy,2);
ywo_e = [ywo_m, ywo_m] + 2*[-std(woy,[],2),std(woy,[],2)];

% Plot the null distribution and the empirical difference
figName = sprintf('Test_PROB_DET_model_rmse_mean_HIST');
fh = mrvNewGraphWin(figName);
patch([xhis,xhis],y_e(:),histcolor{1},'FaceColor',histcolor{1},'EdgeColor',histcolor{1}); % Distribution as the +/- 2SD
hold on
patch([woxhis,woxhis],ywo_e(:),histcolor{2},'FaceColor',histcolor{2},'EdgeColor',histcolor{2}); % Distribution as the +/- 2SD
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',h1.ylim, ... 
        'xlim',h1.xlim, ...
        'ytick',h1.ytick, ...
        'xtick',h1.xtick, ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('rmse','fontsize',16')
% (3) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
dprime_all_voxels = diff([mean(nullDistributionP,1);mean(nullDistributionD,1)]) ...
               ./sqrt(sum([std(nullDistributionP,[],1);std(nullDistributionD,[],1)].^2,1));
title(sprintf('Strength of connection evidence %2.3f',mean(dprime_all_voxels)), ...
    'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),'eps')

% Now let's test whether the prediction of the Probabilistic coectome in
% the voxels where the deterministic did not go is different from zero.
% Which is the prediction fo the deterministic conenctome, the one that has
% no fibers in these voxels.

% Voxels in the Probabilistic connectome that are
% also in the Deterministic connectome.
probOnlyVoxleIdx = ~prob.coordsIdx ;

% Make a new vector for the psig of the Deterministic. THis vector will
% contain zeros at each voxel in Probabilistic not contained in
% Deterministic. It will contain det.psig in the rest of the voxles.
%prob_only_psig = zeros(nbvas{2},numel(find(probOnlyVoxleIdx)));
prob_only_psig = psig{2}(:,probOnlyVoxleIdx);

% Now extract the measured diffusion signal the RMSE. NOTE: 
% The measured signal is the one from the Probabilistic connectome 
% (from the second data set).
dSig         = feGetRep(fe,'dsigdemeanedvox');
dSigProbOnly = dSig(:,probOnlyVoxleIdx);

% Compute the RMSE
rmse_prob_only = sqrt( mean((dSigProbOnly - prob_only_psig).^2,1) );

% THis is the null model, namely that the signal in these voxels where the
% detrministic tractography did not go is zero.
det_only_psig = zeros(size(dSigProbOnly));
rmse_det_only = sqrt( mean((dSigProbOnly - det_only_psig).^2,1) );

% Make a scatter plot.
figNameRmse = sprintf('prob_vs_det_rmse_map_prob_only_voxels');
fhRmseMap = mrvNewGraphWin(figNameRmse);
[ymap,x]  = hist3([rmse_det_only;rmse_prob_only]',{[10:1:70], [10:1:70]});
ymap = ymap./length(prob.Mrmse);
sh = imagesc(flipud(log10(ymap)));
cm = colormap(flipud(hot)); view(0,90);
axis('square')      
set(gca, ...
    'xlim',[1 length(x{1})],...
    'ylim',[1 length(x{1})], ...
    'ytick',[1 (length(x{1})/2) length(x{1})], ...
    'xtick',[1 (length(x{1})/2) length(x{1})], ...
    'yticklabel',[x{1}(end) x{1}(round(end/2)) x{1}(1)], ...
    'xticklabel',[x{1}(1)   x{1}(round(end/2)) x{1}(end)], ...
    'tickdir','out','ticklen',[.025 .05],'box','off', ...
    'fontsize',fontSiz','visible','on')
hold on
plot3([1 length(x{1})],[length(x{1}) 1],[max(ymap(:)) max(ymap(:))],'k-','linewidth',1)
ylabel('Deterministic_{rmse}','fontsize',fontSiz)
xlabel('Probabilistic_{rmse}','fontsize',fontSiz)
cb = colorbar;
tck = get(cb,'ytick');
set(cb,'yTick',[min(tck)  mean(tck) max(tck)], ...
    'yTickLabel',round(1000*10.^[min(tck),...
    mean(tck), ...
    max(tck)])/1000, ...
    'tickdir','out','ticklen',[.025 .05],'box','on', ...
    'fontsize',fontSiz','visible','on')
saveFig(fhRmseMap,fullfile(saveDir,figNameRmse),'eps')

% Make a statisitcal test. To show that the Probabilistic model is better
% than the deterministic model when using all the voxels.
sizeWith    = length(det.rmse_prob);
nullDistributionP = nan(nboots,nmontecarlo);
nullDistributionD = nan(nboots,nmontecarlo);

for inm = 1:nmontecarlo
    parfor ibt = 1:nboots
        nullDistributionP(ibt,inm) = mean(randsample(rmse_prob_only, sizeWith,true));      
        nullDistributionD(ibt,inm) = mean(randsample(rmse_det_only, sizeWith,true));
    end
    
    % Distribution With
    [y(:,inm),xhis] = hist(nullDistributionP(:,inm),linspace(h2.xlim(1),h2.xlim(2),200));
    y(:,inm) = y(:,inm)./sum(y(:,inm));
    
    % Distribution without
    [woy(:,inm),woxhis] = hist(nullDistributionD(:,inm),linspace(h2.xlim(1),h2.xlim(2),200));
    woy(:,inm) = woy(:,inm)./sum(woy(:,inm));
end
y_m = mean(y,2);
y_e = [y_m, y_m] + 2*[-std(y,[],2),std(y,[],2)];

ywo_m = mean(woy,2);
ywo_e = [ywo_m, ywo_m] + 2*[-std(woy,[],2),std(woy,[],2)];

% Plot the null distribution and the empirical difference
figName = sprintf('Test_PROB_DET_model_rmse_in_voxel_of_prob_connectome_only_mean_HIST');
fh = mrvNewGraphWin(figName);
patch([xhis,xhis],y_e(:),histcolor{1},'FaceColor',histcolor{1},'EdgeColor',histcolor{1}); % Distribution as the +/- 2SD
hold on
patch([woxhis,woxhis],ywo_e(:),histcolor{2},'FaceColor',histcolor{2},'EdgeColor',histcolor{2}); % Distribution as the +/- 2SD
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',h2.ylim, ... 
        'xlim',h2.xlim, ...
        'ytick',h2.ytick, ...
        'xtick',h2.xtick, ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('rmse','fontsize',16')
% (3) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
dprime_prob_only = diff([mean(nullDistributionP,1);mean(nullDistributionD,1)]) ...
               ./sqrt(sum([std(nullDistributionP,[],1);std(nullDistributionD,[],1)].^2,1));
title(sprintf('Strength of connection evidence %2.3f',mean(dprime_prob_only)), ...
    'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),'eps')

% Make bar graph showing the Effects. The Probabilistic conenctome is
% better than the deterministic one in all voxels. Moreover it is better
% than a null model in the voxels where the Deterministic did not go. THis
% means that the fiber sgenerated by the probabilistic connectome in those
% voxels generate a fair prediction for the signal in the voxels.
% Plot the null distribution and the empirical difference
figName = sprintf('Test_PROB_DET_model_rmse_mean_BAR');
fh = mrvNewGraphWin(figName);
dprime = [dprime_all_voxels;dprime_prob_only]';
mStrength = mean(dprime,1);
eStrength = [mStrength;mStrength] + 2*[-std(dprime,[],1);std(dprime,[],1)];
hb = bar(mStrength,'facecolor','k');
hold on
plot([1 1; 2 2]',[eStrength],'r-','linewidth',4)
set(gca,'xlim',[0.5 2.5],...
    'ylim',    [0 bymax], ...
    'ytick',[0 bymax/2 bymax], ...
    'xtick',[1 2], ...
    'xticklabel',{'All voxels','only prob voxels'}, ...
    'tickdir','out','box','off', ...
    'fontsize',16,'visible','on');
ylabel('S','FontSize',16)
saveFig(fh,fullfile(saveDir,figName),'eps')

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [algo, lmax, bval, diffModParamsType] = getConditions(runType)
%
% trackingType,lmax,bval,rep,diffusionModelParams
%
% Returns the desired options for a tractography.
switch runType
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 2000
    case 1
        algo = 'p';lmax = 2;bval=2000;diffModParamsType = [1,0];
    case 2
        algo = 'p';lmax = 4;bval=2000;diffModParamsType = [1,0];
    case 3
        algo = 'p';lmax = 6;bval=2000;diffModParamsType = [1,0];
    case 4
        algo = 'p';lmax = 8;bval=2000;diffModParamsType = [1,0];
    case 5
        algo = 'p';lmax = 10;bval=2000;diffModParamsType = [1,0];
    case 6
        algo = 'p';lmax = 12;bval=2000;diffModParamsType = [1,0];
    case 7
        algo = 'p';lmax = 14;bval=2000;diffModParamsType = [1,0];
    case 8
        algo = 'p';lmax = 16;bval=2000;diffModParamsType = [1,0];
      
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 4000
    case 9
        algo = 'p';lmax = 2;bval=4000;diffModParamsType = [1,0];
    case 10
        algo = 'p';lmax = 4;bval=4000;diffModParamsType = [1,0];
    case 11
        algo = 'p';lmax = 6;bval=4000;diffModParamsType = [1,0];
    case 12
        algo = 'p';lmax = 8;bval=4000;diffModParamsType = [1,0];
    case 13
        algo = 'p';lmax = 10;bval=4000;diffModParamsType = [1,0];
    case 14
        algo = 'p';lmax = 12;bval=4000;diffModParamsType = [1,0];
    case 15
        algo = 'p';lmax = 14;bval=4000;diffModParamsType = [1,0];
    case 16
        algo = 'p';lmax = 16;bval=4000;diffModParamsType = [1,0];
         
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 1000
    case 17
        algo = 'p';lmax = 2;bval=1000;diffModParamsType = [1,0];
    case 18
        algo = 'p';lmax = 4;bval=1000;diffModParamsType = [1,0];
    case 19
        algo = 'p';lmax = 6;bval=1000;diffModParamsType = [1,0];
    case 20
        algo = 'p';lmax = 8;bval=1000;diffModParamsType = [1,0];
    case 21
        algo = 'p';lmax = 10;bval=1000;diffModParamsType = [1,0];
    case 22
        algo = 'p';lmax = 12;bval=1000;diffModParamsType = [1,0];
    case 23
        algo = 'p';lmax = 14;bval=1000;diffModParamsType = [1,0];
    case 24
        algo = 'p';lmax = 16;bval=1000;diffModParamsType = [1,0];
         
    % MRTRIX deterministic tractogrpahy 
    % BVAL = 2000
    case 1+16
        algo = 'd';lmax = 2;bval=2000;diffModParamsType = [1,0];
    case 2+16
        algo = 'd';lmax = 4;bval=2000;diffModParamsType = [1,0];
    case 3+16
        algo = 'd';lmax = 6;bval=2000;diffModParamsType = [1,0];
    case 4+16
        algo = 'd';lmax = 8;bval=2000;diffModParamsType = [1,0];
    case 5+16
        algo = 'd';lmax = 10;bval=2000;diffModParamsType = [1,0];
    case 6+16
        algo = 'd';lmax = 12;bval=2000;diffModParamsType = [1,0];
    case 7+16
        algo = 'd';lmax = 14;bval=2000;diffModParamsType = [1,0];
    case 8+16
        algo = 'd';lmax = 16;bval=2000;diffModParamsType = [1,0];
      
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 4000
    case 9+16
        algo = 'd';lmax = 2;bval=4000;diffModParamsType = [1,0];
    case 10+16
        algo = 'd';lmax = 4;bval=4000;diffModParamsType = [1,0];
    case 11+16
        algo = 'd';lmax = 6;bval=4000;diffModParamsType = [1,0];
    case 12+16
        algo = 'd';lmax = 8;bval=4000;diffModParamsType = [1,0];
    case 13+16
        algo = 'd';lmax = 10;bval=4000;diffModParamsType = [1,0];
    case 14+16
        algo = 'd';lmax = 12;bval=4000;diffModParamsType = [1,0];
    case 15+16
        algo = 'd';lmax = 14;bval=4000;diffModParamsType = [1,0];
    case 16+16
        algo = 'd';lmax = 16;bval=4000;diffModParamsType = [1,0];
         
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 1000
    case 17+16
        algo = 'd';lmax = 2;bval=1000;diffModParamsType = [1,0];
    case 18+16
        algo = 'd';lmax = 4;bval=1000;diffModParamsType = [1,0];
    case 19+16
        algo = 'd';lmax = 6;bval=1000;diffModParamsType = [1,0];
    case 20+16
        algo = 'd';lmax = 8;bval=1000;diffModParamsType = [1,0];
    case 21+16
        algo = 'd';lmax = 10;bval=1000;diffModParamsType = [1,0];
    case 22+16
        algo = 'd';lmax = 12;bval=1000;diffModParamsType = [1,0];
    case 23+16
        algo = 'd';lmax = 14;bval=1000;diffModParamsType = [1,0];
    case 24+16
        algo = 'd';lmax = 16;bval=1000;diffModParamsType = [1,0];
        
    % MRTRIX tensor-based tractogrpahy 
    % BVAL = 2000,4000,1000
    case 25+16
        algo = 't';lmax = 2;bval=2000;diffModParamsType = [1,0];
    case 26+16
        algo = 't';lmax = 2;bval=4000;diffModParamsType = [1,0];
    case 27+16
        algo = 't';lmax = 2;bval=1000;diffModParamsType = [1,0];
  
  otherwise
        keyboard
end

end
