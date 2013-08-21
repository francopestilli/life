function fe = s_fe_fit()
%
% This function llustrates how to:
%  - initialize a LIFE structure from a candidate connectome
%  - Generate an optimized connectome from a cadidate LIFE structure
%
%  fe = s_fe_fit()
% 
% INPUTS:  none
% OUTPUTS: fe structure the optimized life structure
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

% Get the base directory for the data
datapath = pestilliDataPath;

dwiFile       = fullfile(datapath,'diffusion','subject1_scan1_2mm_150dir_b2000.nii.gz');
dwiFileRepeat = fullfile(datapath,'diffusion','subject1_scan2_2mm_150dir_b2000.nii.gz');
t1File        = fullfile(datapath,'anatomy','subject1_t1_anatomy.nii.gz');
fgFileName    = fullfile(datapath,'connectomes','subject1_fibers_unculled_2mm_150dir_b2000_probabilistic_lmax8.mat');
savedir       = datapath;
feFileName    = 'fe_culled_subject1_scan1_2mm_150dir_b2000';

% Intialize a local matlab cluster if the parallel toolbox is available.
feOpenLocalCluster;

% Initialize the Connectome
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,savedir,dwiFileRepeat,t1File);

% Now reduce the size of the fiber groups
fe = feConnectomeCull(fe);

% Save it
feConnectomeSave(fe);

return
