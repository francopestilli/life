function opts = feFiberTrackOpts(varargin)
%
% opts = feFiberTrackOpts([options])
%
% TrackOpts is a struct with the following fields:
%  
%  % TEND and FACT parameters
%  opts.stepSizeMm       = 1;
%  opts.faThresh         = 0.20;
%  opts.lengthThreshMm   = 20;
%  opts.angleThresh      = 30; % No turns more than 30 deg in angle
%  opts.wPuncture        = 0.2;
%  opts.whichAlgorithm   = 1;
%  opts.whichInterp      = 1;
%  opts.seedVoxelOffsets = [-0.25  0.25];
%  opts.offsetJitter     = 0;
%
%  % MRTRIX parameters
%  opts.lmax   = 10;
%  opts.nSeeds = 500000;
%
% See Also: dtiFiberTracker, dtiFiberTractOpts .
%
% Franco (c) Stanford Vista Team, 

opts.stepSizeMm       = .2;
opts.faThresh         = 0.15;
opts.lengthThreshMm   = 20;
opts.angleThresh      = 30; % No turns more than 30 deg in angle
opts.wPuncture        = 0.2;
opts.whichAlgorithm   = 1;
opts.whichInterp      = 1;
opts.seedVoxelOffsets = [-0.125 0.125];
opts.offsetJitter     = 0.0;
opts.lmax             = 10;
opts.nSeeds           = 500000;

if(~isempty(varargin))
    for(ii=1:2:numel(varargin)-1)
        opts = setfield(opts, varargin{ii}, varargin{ii+1});
    end
end
        
return;
