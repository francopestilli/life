function fe = feConnectomeSelectFibers(fe, fibersToKeep,recomputeFibersInfo)
% Deletes a set of fibers from the signal and the M matrix.
%
%   fe = feConnectomeSelectFibers(fe,fibersToKeep,[recomputeFibersInfo=0])
%
% The connectome select process builds the parameters inside of fe.life for
% a selected set of fibers.  The parameters inside of fe.life that are
% modified are
%
%   The M matrix
%   The fe.life.fibers entries
%
% We store and keep the whole connectome in .fg.  But not all calculations
% with the M matrix use all the fibers.  In some cases we run smaller
% subsets.
%
% Inputs:
%   fe                  - an fe structure, see feCreate.m, and v_lifeExample.m
%   fibersToKeep        - a list of indexes to the fibers to keep e.g., [1 10 100].
%   recomputeFibersInfo - recomputes the fibers' info this takes a long time. 
%                         So the default is not to recompute (0).
%
% Example:
%   fibersToKeep = 1:50;
%   feConnectomeSelectFibers(fe, fibersToKeep)
%
% See also:  fgExtract
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com.

if notDefined('recomputeFibersInfo'),recomputeFibersInfo=0;end
fprintf('[%s], Selecting fibers from the full Connectome...\n',mfilename)

% Change the name of the fe structure
indx       = strfind(feGet(fe,'name'),'-');
feFileName = feGet(fe,'name');
if ~isempty(indx)
  feFileName = feFileName((indx(end) + 1):end); % Strip away any previous date-string tag.
end
feFileName = sprintf('%s-%s', datestr(now,30),feFileName);
fe         = feSet(fe, 'name',feFileName);

% Delete fibers' columns from the model.
fe = feSet(fe,'Mfiber',feGet(fe,'mkeepfibers',fibersToKeep));

% Remove the fibers from the fiber group:
fe = feSet(fe,'fg',feGet(fe,'fiberssubset',fibersToKeep));

% Take care of the field fe.life.fibers
% Set the subset of tensors, the one only for the left-over fibers
if ~isempty(fe.life.fibers.tensors)
  fe = feSet(fe,'tensors',feGet(fe,'tensors',fibersToKeep));
end

% Remove any fit or cross-validation
fe = feSet(fe,'fit', []);
fe = feSet(fe,'xvalfit',[]);

% By default we do not recompute the fibers info.
if recomputeFibersInfo
  % We disregard fibers that have identical trajectories within the ROI.
  roi = feGet(fe,'roi coords');
  fe  = feSet(fe,'voxel 2 fiber node pairs',fefgGet(feGet(fe,'fg'),'v2fn',roi));
  fe  = feGetConnectomeInfo(fe);
else
  % this is done automatically by feSet: 
  % fe = feSet(fe,'voxel 2 fiber node pairs',[]);
  fe = feSet(fe,'index to unique fibers in each voxel',[]);
  fe = feSet(fe,'number of unique fibers in each voxel', []);
  fe = feSet(fe,'number of total fibers in each voxel', []);
  fe = feSet(fe,'index of total fibers in each voxel', []);
end

return

