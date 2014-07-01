function [fe, indicesFibersKept] = feConnectomeReduceVoxels(fe,voxelsToKeep)
% Select the voxels to keep (extract) in a connectome matrix
%
%   [fe, indicesFibersKept] = feConnectomeReduceVoxels(fe,voxelToKeep)
%
% voxelsToKeep:  Is a binary list of voxels we preserve.
%
% We expand the voxelsToKeep into a binary list of 0's and 1's that in
% which each is expanded by nBvecs.  The 1s are the rows of the M matrix we
% will keep.
%
%  Example:
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com.

% Get the indices to each voxels' signal
vxRows = feGet(fe,'voxelrows',voxelsToKeep);

% Return only the mode and the signal for the voxels we want to keep
fe.life.Mfiber = fe.life.Mfiber(vxRows,:);
fe.life.dSig   = fe.life.dSig(vxRows);

% Set the new number of voxels, by indexing inside the roi and
% returning it as an ROI.
fe.roi.coords = feGet(fe,'roi coords subset',voxelsToKeep);

% Set the new diffusion signal, the one for only these subset of voxels.
fe.life.diffusion_signal_img = fe.life.diffusion_signal_img(voxelsToKeep,:);

% Set the diffusion signal at 0 diffusion weighting (B0) for this voxel:
if ~isempty(fe.life.diffusion_S0_img)
fe.life.diffusion_S0_img = fe.life.diffusion_S0_img(voxelsToKeep);
end

% Set the voxels to fiber/node pairs for a subset of voxels in the conncetome.
fe.life.vovel2FNpair = [];

% Now remove singals for the second data set if it was loaded
if isfield(fe,'rep')   
    if ~isempty(fe.rep.diffusion_signal_img)
    % Set the new diffusion signal, the one for only these subset of voxels.
    fe.rep.diffusion_signal_img = fe.rep.diffusion_signal_img(voxelsToKeep,:);
    end
    
    if ~isempty(fe.rep.diffusion_S0_img)
    % Set the diffusion signal at 0 diffusion weighting (B0) for this voxel:
    fe.rep.diffusion_S0_img = fe.rep.diffusion_S0_img(voxelsToKeep); 
    end
    
end

% Now that we have removed some voxels fom the model, we need to remove also
% the fibers that do not go through the coordinates left in the roi of the model.
% These fibers make no contribution to the signal in the voxels.
%
% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end
% Find the unique fibers in the new ROI.
fibersToKeep = feGet(fe,'uniquefibersindicesinroi');

% Find the indices of the fibers that were deleted
indicesFibersKept = zeros(size(feGet(fe,'fiber weights')));
indicesFibersKept(fibersToKeep) = 1;

if ~poolwasopen, matlabpool close; end

% Remove the fibers of the fascicle from the fe.
fe = feConnectomeReduceFibers(fe, fibersToKeep );


return
