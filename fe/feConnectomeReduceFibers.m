function fe = feConnectomeReduceFibers(fe, fibersToKeep)
% Deletes a set of fibers from the signal and the M matrix.
%
%   fe = feConnectomeDeleteFibers(fe,fibersToKeep)
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
% Franco (c) Stanford Vista Team 2012

if all(fibersToKeep==0)
    fe.life.Mfiber = sparse(size(fe.life.Mfiber,1),size(fe.life.Mfiber,2));
else
    % Delete fibers' columns from the model
    fe.life.Mfiber = fe.life.Mfiber(:,find(fibersToKeep));
end

% Clear the fields that depend o the original fiber group. These are:
% (1) The fit of the model.
if isfield(fe.life,'fit') && ~isempty(fe.life.fit)
    fe.life.fit = [];
end
if isfield(fe.life,'voxfit') && ~isempty(fe.life.voxfit)
    fe.life.voxfit = [];
end
% The field 'fibers' containing some statistics obtained from the original
% fg
if isfield(fe.life,'fibers') && ~isempty(fe.life.fibers.tensors)
    fe.life.fibers.tensors = [];
    fe.life.fibers.total   = [];
    fe.life.fibers.unique  = [];
end

% The actual fiber group.
if all(fibersToKeep==0)
    fe.fg = fgCreate('name',fe.fg.name);
else
    fe.fg = fgExtract(feGet(fe,'fibers img'),find(fibersToKeep),'keep');
end

return

