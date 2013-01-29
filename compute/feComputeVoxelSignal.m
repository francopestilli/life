function voxDSig = feComputeVoxelSignal(fe,voxIndex)
% Predict the signal in a single voxel given the fiber going through the
% voxel.
%
%   voxDSig = feComputeVoxelSignal(fe,voxelIndex)
%
% Franco (c) Stanford Vista Team 2012

% Extract information regarding, voxels, signal and fibers.
S0                = feGet(fe,'b0signalimage',voxIndex);     % non diffusion-weighted signal
bvecs             = feGet(fe,'bvecs');                      % bvecs
bvals             = feGet(fe,'bvals');                      % bvals
tot_fibers_num    = feGet(fe,'tot f num',      voxIndex);   % number of total fibers in the voxel
unique_fibers_num = feGet(fe,'unique f num',     voxIndex); % number of unique fibers in the voxel
tot_fiber_index   = cell2mat(feGet(fe,'totf', voxIndex));   % indexes to the total fibers in the voxels
unique_fiber_index= cell2mat(feGet(fe,'uniquef',voxIndex)); % indexes to the unique fibers in the voxels
voxTensors        = feGet(fe,'voxeltensors',voxIndex);      % Get the tensors for each node in each fiber 
                                                            % going through this voxel

% Compute the predicted signal by each tensors of each node in this voxel.
voxDSig = feComputeSignal(S0, bvecs, bvals, voxTensors);

% Combine the diffusion predictions across nodes of a single fiber.
% Use only the prediction from the unique fibers, not from all the fibers.
if tot_fibers_num ~= unique_fibers_num
  combineM = zeros(tot_fibers_num, unique_fibers_num);
  
  % The matrix combineM is a set of 0s and 1s that will sum together the
  % nodes from a single fiber.
  for ii=1:unique_fibers_num
    combineM(:,ii) = ( tot_fiber_index == unique_fiber_index(ii) );
  end
  
  % The matrix for this voxel starts with each node, and when we multiply
  % by combineM. The resulting matrix represents each fiber
  % (not each node) as a column
  voxDSig = voxDSig*combineM;
end

return
