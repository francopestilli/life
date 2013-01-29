function [fd nodesInVoxel] = mctComputeFiberDensity(fg,roiCoords,varargin)
% 
% function [fd nodesInVoxel] = mctComputeFiberDensity(fg,roiCoords,varargin)
%
% Compute fiber density.
%
% The fiber density is the sum of the weigths per fiber-nodes in eahc
% voxels.
%
% If the weights are all ones (default), the fiber density is the total
% number of nodes in a voxel.
%
% See also: mctBuildDiffusionModel.m, mctFitDiffusionModel.m, fefgGet
%
% Example:
%  See s_mct_fact_cc_prefrontal_rois.m
%   [A dSig] = mctDiffusionModel(dwi,fgImg,theseCoords);
%    w       = mctDiffusionModeFit(A,dSig);
%    fdBeforeMct = mctComputeFiberDensity(fg,roiCoords);
%    fdAfterMct = mctComputeFiberDensity(fg,roiCoords,w);
%
% (c) 2012 Stanford VISTA Team


% se tdefault for optional parameters
voxel2FNpair = []; weights = [];
optargin = size(varargin,2);
if (optargin >= 1), weights = varargin{1};end
if (optargin >= 2), voxel2FNpair = varargin{2};end

% if the weights were not passed in we assume that there are equal
% weights across fibers
if isempty(weights)
    weights = ones(size(fg.fibers));
end

% if the voxel2FNpar cell was not passed in we compute it
% this is slow though...
if isempty(voxel2FNpair)
    voxel2FNpair = fefgGet(fg,'voxel 2 fiber node pairs',roiCoords);
end

% compute the fiber density for each voxel
which_fibers = weights > 0;
for ii = 1:length(voxel2FNpair) % number of voxels
    % compute how many nodes are in each voxel
    nodesInVoxel(ii) = sum( which_fibers( voxel2FNpair{ii}(:,1) ) );
    
    % compute the expected fascicle density 
    fd(ii) = sum( weights( voxel2FNpair{ii}(:,1) ) );
end

return