function s_ms_batch_multi_track_proclus_SUBJID(runType)
%
% Tracks using MRTRIX with different parameters.
%
% INPUTS
%  runType - is a number that indexes into the specific settings saved
%  inside the local function getContions.
%
% Franco (c) Stanford Vista Team, 2013

% Select the parameters for the current conditions
[trackingAlgorithm, dtFile,fibersFolder,lmax,nSeeds,wmMask] = getConditions(runType);

% Run the tracking condition
[status, results, fg, pathstr] = feTrack(trackingAlgorithm, dtFile,fibersFolder,lmax,nSeeds,wmMask);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [algo, dtFile,fibersFolder,lmax,nSeeds,wmMask] = getConditions(runType)

algo         = {'probabilistic'};
nSeeds       = 500000;
dtFile       = '/biac2/wandell2/data/diffusion/PATH/to/dt6.mat';
wmMask       = [];
lmax         = 2:2:12;

switch runType
  %% Subject 1
  case 1
    fibersFolder = '/azure/scr1/htakemura/connectomes/SUBJ1';
    
    %% subject 2
  case 2
    fibersFolder = '/azure/scr1/htakemura/connectomes/SUBJ2'; 
    
  otherwise
    keyboard
end
end
