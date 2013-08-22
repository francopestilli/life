function fe = s_fe_fit()
%
% This function llustrates how to:
%  - initialize a LIFE structure from a candidate connectome
%  - Generate an optimized connectome from a cadidate connectome using the
%  LIFE strustrue
%
%  fe = s_fe_fit()
% 
% INPUTS:  none
% OUTPUTS: fe structure the optimized life structure
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

% Get the base directory for the data
datapath = pestilliDataPath;

% Build the file names for the diffusion data, the anatomical MR, the fiber
% group containing the connectome and the 
dwiFile       = fullfile(datapath,'diffusion','subject1_scan1_2mm_150dir_b2000.nii.gz');
dwiFileRepeat = fullfile(datapath,'diffusion','subject1_scan2_2mm_150dir_b2000.nii.gz');
t1File        = fullfile(datapath,'anatomy','subject1_t1_anatomy.nii.gz');
fgFileName    = fullfile(datapath,'connectomes','subject1_fibers_unculled_2mm_150dir_b2000_probabilistic_lmax8.mat');
savedir       = datapath;

% The final connectome and dat astructure will be saved with this name:
feFileName    = 'fe_culled_subject1_scan1_2mm_150dir_b2000';

% Intialize a local matlab cluster if the parallel toolbox is available.
feOpenLocalCluster;

% Initialize the Connectome
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,savedir,dwiFileRepeat,t1File);

% Fit the model and cull. This will take some time...
fe = feConnectomeCull(fe);

% Fastest way to cull the connectome
% fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

% Save it
feConnectomeSave(fe);

% Make a plot of the weights:
w = feGet(fe,'fiber weights');
figName = sprintf('Fascicle weights');
mrvNewGraphWin(figName);
[y,x] = hist(w(w>0),logspace(-4,-.3,50));
semilogx(x,y)
set(gca,'tickdir','out','fontsize',16,'box','off')
title('fascicle weights','fontsize',16)
ylabel('number of fascicles','fontsize',16)
xlabel('fascicle weight','fontsize',16)

% Make a plot of the RMSE:
rmse   = feGet(fe,'vox rmse');
rmsexv = feGetRep(fe,'vox rmse');
figName = sprintf('RMSE');
mrvNewGraphWin(figName);
% Non-cross-validated
[y,x] = hist(rmse,50);
plot(x,y,'k-')
hold on
% Cross-validated
[y,x] = hist(rmsexv,50);
plot(x,y,'r-')
set(gca,'tickdir','out','fontsize',16,'box','off')
title('Root-mean squared error distribution across voxels','fontsize',16)
ylabel('number of voxels','fontsize',16)
xlabel('rmse','fontsize',16)
legend({'RMSE fitted data set','RMSE cross-validated'},'fontsize',16)

% Make a plot of the RMSE Ratio:
R   = feGetRep(fe,'voxrmseratio');
figName = sprintf('RMSE RATIO');
mrvNewGraphWin(figName);
[y,x] = hist(R,linspace(.5,2,50));
plot(x,y,'k-')
hold on
plot([1 1],[0 450],'k--')
set(gca,'tickdir','out','fontsize',16,'box','off')
title('Root-mean squared error Ratio','fontsize',16)
ylabel('number of voxels','fontsize',16)
xlabel('R_{rmse}','fontsize',16)


return

