function [status, results, fg, pathstr] = feTrack(trackingAlgorithm, dtFile,fibersFolder,lmax,nSeeds,wmMask)
%
% This function creates whole-brain white-matter connectome.
%
% Inputs:
%  - dt6Dir - directory containing the dt6.mt file.
%  - fibersFolder - folder to use to save all the fibers computed.
%  - lmax - Max harmonic order. Default = 10.
%  - nSeeds - number of seeds to use for mrtrix, this determines the final
%             number of fibers.
%
% Franco (c) Stanford VISTA Team 2012

if notDefined('trackingAlgorithm'), trackingAlgorithm = {'tensor'};end
if notDefined('dtFile'),            dtFile = ...
    '/biac2/wandell6/data/frk/LiFE/data/fp20120420/150dirs_b1000_1/dt6.mat';end
if notDefined('wmMask'),            wmMask  = [];end
if notDefined('fibersFolder'),      fibersFolder = ...
    '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep1';end
if notDefined('lmax'),              lmax    = [8];end
if notDefined('nSeeds'),            nSeeds  = 500000; end
  
runInBackground = 0;
verbose = 1;
  
% make the fibers' folder if it does not exist
if ~exist(fibersFolder,'dir'), system(sprintf('mkdir -v %s',fibersFolder)); end

fprintf('\n\n[%s] Performing whole-brain tractograpy ...\n\n',mfilename);

% Track with MRTRIX
for il = 1:length(lmax)
  % This first step initializes all the files necessary for mrtrix.
  % This can take a long time.
  files           = mrtrix_init(dtFile, lmax(il),fibersFolder,wmMask);

  
  % MRTRIX - We run this first because it is fast.
  for ii = 1:length(trackingAlgorithm)
    % Track and save fibers using mrtrix
    [status, results, fg, pathstr] = mrtrix_track(files, files.brainmask, files.wm, switchAlgo(trackingAlgorithm{ii}), nSeeds, runInBackground, verbose);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%
% Choose algorithm. %
%%%%%%%%%%%%%%%%%%%%%
function algo = switchAlgo(algo)
% Gets the algorithm type for tracktography.
%
% Tracking algorithm: 
% 
% 1=prob, 2=stream, 3=tensor.
%
% Franco

if isnumeric(algo) % from index to string
  switch algo
    case {1}
      algo = 'prob';
    case {2}
      algo = 'stream';
    case {3}
      algo = 'tensor';
    otherwise
      keyboard
  end
  
elseif ischar(algo)
  switch algo
     case {'prob','probabilistic','mrtrix probabilistic'}
      algo = 'prob';
    case {'stream','streamline','deterministic', 'mrtrix deterministic'}
      algo = 'stream';
    case {'tensor','tensor based','mrtrix tensor based tractography'}
      algo = 'tensor';
    otherwise
      keyboard
  end
  
else
  keyboard
end
