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

switch runType
  %% Subject 1
  case 1
    fibersFolder = '/azure/scr1/frk/FP_96dirs_b2000_1p5iso/life_mrtrix_rep1/';
    algo         = {'probabilistic','deterministic'};
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120718_2975/96dirs_b2000_1point5iso_1/dt6.mat';
    wmMask       = [];
    lmax         = 2:2:12;
  case 2
    fibersFolder = '/azure/scr1/frk/FP_96dirs_b2000_1p5iso/life_mrtrix_rep2/';
    algo         = {'probabilistic','deterministic'};
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120718_2975/96dirs_b2000_1point5iso_1/dt6.mat';
    wmMask       = [];
    lmax         = 2:2:12;
  case 3
    fibersFolder = '/azure/scr1/frk/FP_96dirs_b2000_1p5iso/life_mrtrix_rep3/';
    algo         = {'probabilistic','deterministic'};
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120718_2975/96dirs_b2000_1point5iso_1/dt6.mat';
    wmMask       = [];
    lmax         = 2:2:12;
    
    %% tensor
  case 4
    fibersFolder = '/azure/scr1/frk/FP_96dirs_b2000_1p5iso/life_mrtrix_rep1/';
    algo         = {'tensor'};
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120718_2975/96dirs_b2000_1point5iso_1/dt6.mat';
    wmMask       = [];
    lmax         = 1;
    
  case 5
    fibersFolder = '/azure/scr1/frk/FP_96dirs_b2000_1p5iso/life_mrtrix_rep2/';
    algo         = {'tensor'};
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120718_2975/96dirs_b2000_1point5iso_1/dt6.mat';
    wmMask       = [];
    lmax         = 1;
    
  case 6
    fibersFolder = '/azure/scr1/frk/FP_96dirs_b2000_1p5iso/life_mrtrix_rep3/';
    algo         = {'tensor'};
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120718_2975/96dirs_b2000_1point5iso_1/dt6.mat';
    wmMask       = [];
    lmax         = 1;
   
  otherwise
    keyboard
end
end
