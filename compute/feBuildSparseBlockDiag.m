function A = feBuildSparseBlockDiag(nBvecs,nVoxels, values)
%  Build a sparse block-diagonal matrix with ones in all filled elements.
%
%       A = feBuildSparseBlockDiag(nBvecs,nVoxels, [values]);
%
%  This is used to build the matrix used to represent the isotropic
%  component contributions to the diffusion data in each voxel.
%
%  INPUTS:
%  nBvecs  -  How many directions in your DWI data.
%  nVoxels -  How many voxels in the ROI.
%  values  -  Either a value for every voxel (column, size(values) = nVoxels) OR
%             one value for all the voxels, e.g., 0.5
%
%  OUTPUTS:
%  A       - A block-diagonal matrix.
%            In each column, there are nBvecs filled rows. These are always filled
%            with either a '1' OR a value (if passed in). The location of the
%            non-zero elements in each column is:
%            col * nBvecs + 1: (col + 1) * nBvecs
%
%  EXAMPLE:
%    >> A = feBuildSparseBlockDiag(3,3)
%       A = (1,1) 1
%           (2,1) 1
%           (3,1) 1
%           (4,2) 1
%           (5,2) 1
%           (6,2) 1
%           (7,3) 1
%           (8,3) 1
%           (9,3) 1
%
%  (c) Stanford VISTA Team, 2012

if notDefined('values'),values = ones(1,nVoxels);end

% If only one value is passed in (the same value for all voxels)
% We replicate it.
if (length(values) == 1), values = values .* ones(1,nVoxels);end

idx_rows    = 1:nBvecs * nVoxels;
idx_columns = reshape(repmat(1:nVoxels,nBvecs,1), 1, nBvecs * nVoxels);
A           = sparse(idx_rows, idx_columns, ones(nBvecs, nVoxels) .* ...
                                            repmat(values,nBvecs,1));

end