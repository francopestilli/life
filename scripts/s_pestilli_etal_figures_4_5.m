function [fh, fe] = s_pestilli_etal_Figures_4_5()
%% Example of initialization and fitting of the LiFE model
%
% This function illustrates how:
%  - Set up an fe (fascicle evaluation) structure.
%  - Combines them with fibers tracted between the two ROIs 
%  - Generates an optimized connectome from a candidate connectome using LIFE 
%  - Performs a virtual lesion
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com.

% Intialize a local matlab cluster if the parallel toolbox is available.
feOpenLocalCluster;

%% Build the file names for the diffusion data, the anatomical MR.
dwiFile       = fullfile(lifeDemoDataPath('diffusion'),'pestilli_etal_life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz');
dwiFileRepeat = fullfile(lifeDemoDataPath('diffusion'),'pestilli_etal_life_demo_scan2_subject1_b2000_150dirs_stanford.nii.gz');
t1File        = fullfile(lifeDemoDataPath('anatomy'),  'pestilli_etal_life_demo_anatomy_t1w_stanford.nii.gz');

%% (1) Probabilistic CSD-based connectome:
% We will analyze first the CSD-based probabilistic tractography
% connectome.
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
                'pestilli_et_al_life_demo_mrtrix_csd_lmax10_probabilistic.mat');

% The final connectome and data astructure will be saved with this name:
feFileName    = 'life_build_model_demo_CSD_PROB';


% Initialize the Connectome.
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFileRepeat,t1File);

% Fit the model.
fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

% Plot a histogram of the RMSE.
[fh(1), ~, ~] = plotHistRMSE(fe,'Probabilistic');

% Plot a histogram of the RMSE ratio.
% The Rrmse is the ratio between data test-retest reliability and 
% model error (the quality of the model fit).
[fh(2), ~] = plotHistRrmse(fe,'Probabilistic');

% Plot a histogram of the fitted fascicle weights. 
[fh(3), ~] = plotHistWeigths(fe,'Probabilistic');

%% (2) Deterministic tensor-based connectome:
% We will analyze first the CSD-based probabilistic tractography
% connectome.
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
                'pestilli_et_al_life_demo_mrtrix_tensor_deterministic.mat');

% The final connectome and data astructure will be saved with this name:
feFileName    = 'life_build_model_demo_TENSOR_DET';

% Initialize the Connectome.
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFileRepeat,t1File);

% Fit the model.
fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

% Plot a histogram of the RMSE.
[fh(1), ~, ~] = plotHistRMSE(fe,'Deterministic');

% Plot a histogram of the RMSE ratio.
% The Rrmse is the ratio between data test-retest reliability and 
% model error (the quality of the model fit).
[fh(2), ~] = plotHistRrmse(fe,'Deterministic');

% Plot a histogram of the fitted fascicle weights. 
[fh(3), ~] = plotHistWeigths(fe,'Deterministic');

keyboard
end

% ---------- Local  Functions ----------- %
function [fh, rmse, rmsexv] = plotHistRMSE(fe,tractograpy)
% Make a plot of the RMSE:

% Extract the RMSE of the model on the fitted data set
rmse   = feGet(fe,'vox rmse');

% Extract the RMSE of the model on the second data set (cross-validated
% error)
rmsexv = feGetRep(fe,'vox rmse');

figName = sprintf('%s - RMSE',tractograpy);
fh = mrvNewGraphWin(figName);
[y,x] = hist(rmse,50);
plot(x,y,'k-');
hold on
[y,x] = hist(rmsexv,50);
plot(x,y,'r-');
set(gca,'tickdir','out','fontsize',16,'box','off');
title('Root-mean squared error distribution across voxels','fontsize',16);
ylabel('number of voxels','fontsize',16);
xlabel('rmse (scanner units)','fontsize',16);
legend({'RMSE fitted data set','RMSE cross-validated'},'fontsize',16);
end

function [fh, R] = plotHistRrmse(fe,tractograpy)
% Make a plot of the RMSE Ratio:
R       = feGetRep(fe,'voxrmseratio');
figName = sprintf('%s - RMSE RATIO',tractograpy);
fh      = mrvNewGraphWin(figName);
[y,x]   = hist(R,linspace(.5,4,50));
plot(x,y,'k-','linewidth',2);
hold on
plot([median(R) median(R)],[0 1200],'r-','linewidth',2);
plot([1 1],[0 1200],'k-');
set(gca,'tickdir','out','fontsize',16,'box','off');
title('Root-mean squared error ratio','fontsize',16);
ylabel('number of voxels','fontsize',16);
xlabel('R_{rmse}','fontsize',16);
legend({sprintf('Distribution of R_{rmse}'),sprintf('Median R_{rmse}')});
end

function [fh, w] = plotHistWeigths(fe,tractograpy)
% Make a plot of the weights:
w       = feGet(fe,'fiber weights');
figName = sprintf('%s - Distribution of fascicle weights',tractograpy);
fh      = mrvNewGraphWin(figName);
[y,x]   = hist(w( w > 0 ),logspace(-5,-.3,40));
semilogx(x,y,'k-','linewidth',2)
set(gca,'tickdir','out','fontsize',16,'box','off')
title( ...
    sprintf('Number of fascicles candidate connectome: %2.0f\nNumber of fascicles in optimized connetome: %2.0f' ...
    ,length(w),sum(w > 0)),'fontsize',16)
ylabel('Number of fascicles','fontsize',16)
xlabel('Fascicle weight','fontsize',16)
end

