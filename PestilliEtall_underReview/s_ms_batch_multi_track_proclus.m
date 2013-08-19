function s_ms_batch_multi_track_proclus(runType)
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
  case 1
    algo         = {'tensor'};
    fibersFolder = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep1/';
    lmax         = 'na';
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120420_2290/150dirs_b1000_1/dt6.mat';
    wmMask       = [];
  case 2
    algo         = {'tensor'};
    fibersFolder = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep1/';
    lmax         = 'na';
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20110922_1125/150dirs_b2000_1/dt6.mat';
    wmMask       = [];
  case 3
    algo         = {'tensor'};
    fibersFolder = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep1/';
    lmax         = 'na';
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120420_2290/150dirs_b4000_1/dt6.mat';
    wmMask       = [];
  case 4
    algo         = {'tensor'};
    fibersFolder = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep2/';
    lmax         = 'na';
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120420_2290/150dirs_b1000_1/dt6.mat';
    wmMask       = [];
  case 5
    algo         = {'tensor'};
    fibersFolder = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep2/';
    lmax         = 'na';
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20110922_1125/150dirs_b2000_1/dt6.mat';
    wmMask       = [];
  case 6
    algo         = {'tensor'};
    fibersFolder = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep2/';
    lmax         = 'na';
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120420_2290/150dirs_b4000_1/dt6.mat';
    wmMask       = [];
  case 7
    algo         = {'tensor'};
    fibersFolder = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep3/';
    lmax         = 'na';
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120420_2290/150dirs_b1000_1/dt6.mat';
    wmMask       = [];
  case 8
    algo         = {'tensor'};
    fibersFolder = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep3/';
    lmax         = 'na';
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20110922_1125/150dirs_b2000_1/dt6.mat';
    wmMask       = [];
  case 9
    algo         = {'tensor'};
    fibersFolder = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep3/';
    lmax         = 'na';
    nSeeds       = 500000;
    dtFile       = '/biac2/wandell2/data/diffusion/pestilli/20120420_2290/150dirs_b4000_1/dt6.mat';
    wmMask       = [];
    
  otherwise
    keyboard
end
end
