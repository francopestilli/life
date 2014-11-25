function [fe, removedFibers] = feConnectomeReduceFibers(fe, fibersToKeep)
% Deletes a set of fibers from the signal and the M matrix.
%
%   [fe, removedFibers] = feConnectomeDeleteFibers(fe,fibersToKeep)
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
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com.

% Return the indices to the fibers removed from the connectome
removedFibers = find(~fibersToKeep);

if all(fibersToKeep==0)
    fe.Mfiber = sparse(size(fe.Mfiber,1),size(fe.Mfiber,2));
else
    % Delete fibers' columns from the model
    fe.Mfiber = fe.Mfiber(:,fibersToKeep);
end

% Clear the fields that depend o the original fiber group. These are:
% (1) The fit of the model.
if isfield(fe.life,'fit') && ~isempty(fe.fit)
    fe.fit.weights         = fe.fit.weights(fibersToKeep);
    fe.fit.results.nParams = sum(fibersToKeep);
end
if isfield(fe.life,'voxfit') && ~isempty(fe.voxfit)
    fe.voxfit = [];
end
% The field 'fibers' containing some statistics obtained from the original
% fg
if isfield(fe.life,'fibers') && ~isempty(fe.fibers.tensors)
    fe.fibers.tensors = [];
    fe.fibers.total   = [];
    fe.fibers.unique  = [];
end

% The actual fiber group.
if all(fibersToKeep==0)
    fe = feSet(fe,'fg img',fgCreate('all fibers were deleted',fe.fg.name,'img'));
else
    % feExtract requires indices to each fibers and does not accept logical
    % inuts. Here we make sure we are passing indices not a logical vector.
    if (length(unique(fibersToKeep)) == 2)
        if   (all(unique(fibersToKeep) == [0, 1]'))
            fibersToKeep = find(fibersToKeep);end
    end
    fe = feSet(fe,'fg img',fgExtract(feGet(fe,'fibers img'),fibersToKeep,'keep'));
end

return

