function feTrack(dt6Dir,fibersFolder,lmax,nSeeds,wmMask)
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

if notDefined('lmax'),         lmax = [22 24];end
if notDefined('dt6Dir'),       dt6Dir = '/biac2/wandell6/data/frk/LiFE/data/fp20120420/150dirs_b1000_1';end
if notDefined('wmMask'),       wmMask = [];end
if notDefined('fibersFolder'), fibersFolder = 'fe_fibers';end
if notDefined('nSeeds'),       nSeeds = 200000; end

% Make the full path to the dt6 for and the folder for the fibers
dtFile  = fullfile(dt6Dir,'dt6.mat');
fibDir  = fullfile(dt6Dir,fibersFolder);

% make the fibers' folder if it does not exist
if ~exist(fibDir,'dir'), system(sprintf('mkdir -v %s',fibDir)); end

fprintf('\n\n[%s] Performing whole-brain tractograpy ...\n\n',mfilename);

% Track with MRTRIX
for il = 1:length(lmax)
  % This first step initializes all the files necessary for mrtrix.
  % This can take a long time.
  files           = mrtrix_init(dtFile, lmax(il),fibersFolder,wmMask);
  runInBackground = 0;
  verbose = 1;
  
  % MRTRIX - We run this first because it is fast.
  for ii = [4,5]
    % Track and save fibers using mrtrix
    % [status, results, fg, pathstr] = mrtrix_track(csd, roi, mask, mode, nSeeds, bkgrnd, verbose)
    mrtrix_track(files.csd, files.brainmask, files.wm, switchAlgo(ii), nSeeds, runInBackground, verbose);
  end
end

% Track with TEND/FACT
%
% Parameters...
% for performing whole-brain tractography using TEND and FACT
opts = feFiberTrackOpts;
for ii = [1,3]
  % set up file name and options
  fgWBname = [];
  fgWBname = fullfile(fibDir,sprintf('fg_whole_brain_%s.mat',switchAlgo(ii)));
  opts.whichAlgorithm = ii;
  
  % track and save the fibers
  feFiberTrackWholeBrain(dtFile,opts,fgWBname,'both',1,'pdb');
end


%%%%%%%%%%%%%%%%%%%%%
% Choose algorithm. %
%%%%%%%%%%%%%%%%%%%%%
function algo = switchAlgo(algo)
% Gets the algorithm type for tracktography given the index stored in
% dtiFiberTraker.m
%
% Tracking algorithm: 
% 0=FACT Euler, 1=FACT RK4, 2=TEND Euler, 3=TEND RK4,
% 4=prob, 5=stream. (Last two use mrtrix)
%
% Franco

if isnumeric(algo) % from index to string  
  switch algo
    case {0}
      algo = 'FACT_euler';
    case {1}
      algo = 'FACT_RK4';
    case {2}  
      algo = 'TEND_euler';
    case {3} 
      algo = 'TEND_RK4';
    case {4}  
      algo = 'prob';
    case {5} 
      algo = 'stream';
    otherwise
      keyboard
  end
  
elseif ischar(index)
  switch algo
    case {'FACTe','FACTE','FACT Euler','FACT_euler'}
      algo = 0;
    case {'FACT','fact rk4','FACT_RK4','FACT RK4'}
      algo = 1;
    case {'TENDe','TENDE','TEND_euler'}
      algo = 2;
    case {'TEND','TEND_RK4','TEND RK4', 'tend rk4'}
      algo = 3;    
    case {'prob','probabilistic','mrtrix probabilistic'}
      algo = 4;
    case {'stream','streamline','determinstic', 'mrtrix deterministic'}
      algo = 5;

    otherwise
      keyboard
  end
  
else
  keyboard
end
