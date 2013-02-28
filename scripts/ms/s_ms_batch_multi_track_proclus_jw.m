function s_ms_batch_multi_track_proclus_jw(runType)
%
% Tracks with different mrtrix parameters.
%
%
%
% Franco (c) Stanford Vista Tema, 2013

% Select the parameters for the current conditions
[trackingAlgorithm, dtFile,fibersFolder,lmax,nSeeds,wmMask] = getConditions(runType);

% Run the tracking condition
[status, results, fg, pathstr] = feTrack(trackingAlgorithm, dtFile,fibersFolder,lmax,nSeeds,wmMask);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [algo, dtFile,fibersFolder,lmax,nSeeds,wmMask] = getConditions(runType)

switch runType
  %% Tensor based
  case 1
    algo         = {'tensor'};
    fibersFolder = '/azure/scr1/frk/JW_96dirs_b2000_1p5iso/life_mrtrix_rep1/';
    lmax         = 'n';
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/winawer/20120410_2202/96dirs_b2000_1point5iso_1/dt6.mat';
    wmMask       = [];
  
  
    %% Probabilistic
  case 2
    algo         = {'probabilistic','deterministic'};
    fibersFolder = '/azure/scr1/frk/JW_96dirs_b2000_1p5iso/life_mrtrix_rep1/';
    lmax         = [2:2:12];
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/winawer/20120410_2202/96dirs_b2000_1point5iso_1/dt6.mat';
    wmMask       = [];
    
 
  otherwise
    keyboard
end
end
