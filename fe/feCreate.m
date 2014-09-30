function fe = feCreate
% Create a linear fascicle evaluation structure
%
%  fe = feCreate
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com.

% Book-keeping
fe.name   = 'default'; % This one's name
fe.type   = 'faseval'; % Always fascicle evaluation

% Overview
fe.fg      = [];  % Fiber group candidate fascicles, uses fgSet/Get
fe.roi     = [];  % Cell array of regions of interest where we evaluate

% Path to files
fe.path.folder = fullfile(pwd,'fe'); % This is the default fascicle evaluation folder location
fe.path.dwifile = [];  % Diffusion weighted file used for testing results
fe.path.dtfile  = [];  % Diffusion weighted file used for testing results
fe.path.savedir = [];  % Top directory under which all files will be saved
fe.path.anatomy = [];  % 3D high-resolution anatomical file
fe.path.dwifilerep = [];  % File used for cross-validating the results. 
                          % A secodn measuremnt ideally acquried in the 
                          % same session witht the same scanning
                          % parameters.

% Life Algorithm parameters
fe.Mfiber       = []; % Fiber portion of A matrix
fe.voxel2FNpair = []; % Voxels to fiber/node pair, fefgGet(fg,'voxel 2 fiber node pair')
fe.dSig         = []; % Signal of fibers alone (demeaned).
fe.fibers       = [];
fe.diffusion_signal_img = [];
fe.diffusion_S0_img     = [];
fe.bvecs                = [];
fe.bvals                = [];
fe.bvecsindices         = [];
fe.imagedim             = [];
fe.xform          = [];  % Transforms between coords for fg, roi, and dwi data
fe.xform.img2acpc = [];
fe.xform.acpc2img = [];

% Repetition file for cross-validation
fe.rep = [];

return
