function fe = feConnectomeReduceVoxels(fe,voxelsToKeep)
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

% GEt the indices to each voxels' signal
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
fe.life.diffusion_S0_img = fe.life.diffusion_S0_img(voxelsToKeep);

% Set the voxels to fiber/node pairs for a subset of voxels in the conncetome.
fe.life.vovel2FNpair = [];

% Now remove singals for the second data set if it was loaded
if isfield(fe,'rep')
    % Set the new diffusion signal, the one for only these subset of voxels.
    fe.rep.diffusion_signal_img = fe.rep.diffusion_signal_img(voxelsToKeep,:);
    
    % Set the diffusion signal at 0 diffusion weighting (B0) for this voxel:
    fe.rep.diffusion_S0_img = fe.rep.diffusion_S0_img(voxelsToKeep); 
end

return