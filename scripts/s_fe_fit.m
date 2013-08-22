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

% Save it
feConnectomeSave(fe);

% Make a plot of the weights:
w = feGet(fe,'fiber weights');
figure
[y,x] = hist(w(w>0),logspace(-8,-.3,50));
semilogx(x,y)
title('fascicle weights')
ylabel('number of fascicles')
xlabel('weight')

% Make a plot of the RMSE:
rmse   = feGet(fe,'vox rmse');
rmsexv = feGetRep(fe,'vox rmse');
figure
% Non-cross-validated
[y,x] = hist(rmse,50);
plot(x,y,'k-')
hold on
% Cross-validated
[y,x] = hist(rmsexv,50);
plot(x,y,'r-')
title('Root-mean squared error distribution across voxels')
ylabel('number of voxels')
xlabel('rmse')
legend({'RMSE fitted data set','RMSE cross-validated'})

% Make a plot of the RMSE Ratio:
R   = feGetRep(fe,'voxrmseratio');
figure
[y,x] = hist(R,linspace(.5,2,50));
plot(x,y,'k-')
hold on
plot([1 1],[0 450],'k--')
title('Root-mean squared error Ratio')
ylabel('number of voxels')
xlabel('R_{rmse}')


return

