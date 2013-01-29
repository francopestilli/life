function [Mfiber dSig] = feConnectomeReduceVoxels(fe,voxelsToKeep)
% Select the voxels to keep (extract) in a connectome matrix
%
%   fe = feConnectomeReduceVoxels(fe,voxelToKeep)
%
% voxelsToKeep:  Is a binary list of voxels we preserve.
%
% We expand the voxelsToKeep into a binary list of 0's and 1's that in
% which each is expanded by nBvecs.  The 1s are the rows of the M matrix we
% will keep.
%
%  Example:
%    Run v_lifeExample to get a fe that is initialized.
%    n = feGet(fe,'n voxels');
%    tmp = sort(unique(randi(n,500,1)));
%    voxelsToKeep = zeros(n,1); voxelsToKeep(tmp) = 1;
%    fe = feConnectomeReduceVoxels(fe,voxelsToKeep);
%
% Franco (c) Stanford Vista Team 2012

Mfiber = fe.life.Mfiber(feGet(fe,'voxelrows',voxelsToKeep),:);
dSig   = fe.life.dSig(feGet(fe,'voxelrows',voxelsToKeep))';

% If we want to actually reduce the fe to a subset of voxels, we would have
% to update the related information inside fe.

% Set the new number of voxels, by indexing inside the roi and
% returning it as an ROI.
%fe.roi.coords = feGet(fe,'roi coords subset',voxelsToKeep);

% Set the new diffusion signal, the one for only these subset of voxels.
%fe.life.diffusion_signal_img = feGet(fe,'dsig invox',voxelsToKeep);

% Set the diffusion signal at 0 diffusion weighting (B0) for this voxel:
%fe.life.diffusion_S0_img = feGet(fe,'s0 vox',voxelsToKeep);

% Set the voxels to fiber/node pairs for a subset of voxels in the conncetome.
%fe.life.vovel2FNpair = feGet(fe,'v2fn sub',voxelsToKeep);

% Take care of the fibers info:
%fe = feGetConnectomeInfo(fe);

return