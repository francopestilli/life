function fe = feConnectomeReduceFibers(fe, fibersToKeep,reduceVoxels)
% Deletes a set of fibers from the signal and the M matrix.
%
%   fe = feConnectomeDeleteFibers(fe,fibersToKeep,[reduceVoxels])
%
% Inputs:
%   - fe, an fe structure, see feCreate.m, and v_lifeExample.m
%   - fibersToKeep, a list of indexes to the fibers to keep e.g., [1 10 100].
%   - reduceVoxels, (optional), if set to 1 it will remove all the voxels
%       from the fe structure where the fibers left inside the fe structure do
%       not have go through.
%
% Example:
%   fibersToKeep = 1:50;
%   feConnectomeReduceFibers(fe, fibersToKeep)
% Franco (c) Stanford Vista Team 2012

if notDefined('reduceVoxels'), reduceVoxels=0;end

% Delete fibers' columns from the model
fe.life.Mfiber = fe.life.Mfiber(:,fibersToKeep);

% Collect the voxels that need to be deleted. By looking at the voxels that
% after removing the fibers have no prediction left in the rows of Mfiber.
usedVoxels = feGet(fe,'used voxels');
nUsed = feGet(fe,'n used voxels');
voxelsToKeep = false(nUsed,1);   % Indicator variable
parfor vv = 1:nUsed   % For each voxel
    
  % Find the signals for every direction (rows) and every kept fiber (cols)
  % in this voxel.
  thisVoxSig = fe.life.Mfiber(feGet(fe,'voxel rows',usedVoxels(vv)),:);  
  
  % If the entries are not all zero, keep the voxel
  if ~isequal(thisVoxSig(:),zeros(numel(thisVoxSig),1))
      voxelsToKeep(vv) = true;
  end
end

% Set the subset of tensors, the one only for the left-over fibers
if ~isempty(fe.life.tensors)
  fe.life.fibers.tensors = feGet(fe,'tensors subset',fibersToKeep);
end

% Take care of the field fe.life.fibers

% Re-set all the fields in fe given the changes in the voxels number given
% by the changes in the fibers number
if reduceVoxels
  fe = feConnectomeReduceVoxels(fe,voxelsToKeep);
end

return

