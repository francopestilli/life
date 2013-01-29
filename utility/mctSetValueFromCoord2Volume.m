function vol = mctSetValueFromCoord2Volume(vals,coords,vol)
%
% function vol = mctSetValueFromCoord2Volume(vals,coords,vol)
%
% Set the values inside a volume specified by some coordinates to the given
% values.
%
% In other terms, this function is meant to be used to change specific
% values inside a volume referred by a set of coordinates.
%
% Inputs:
%    - vals, is a nxm vector of values, one value per voxel in coords.
%    Where n is the number of voxels to be addressed in the volume. 
%    m is the number of measurements (e.g., number of diffusion directions). 
%    - coords, is  nx3 array of image coordinates to
%    index inside vol - vol, is a 3D or 4D volume.
%
% Outputs:
%    - vol, the output has the same size of the iputs vol. All values
%    inside vol are the same except for the one specified in coords, which
%    are set to vals.
%
% Example:
%     vol = mctSetValueFromCoord2Volume(0.5,[1,2,3],ones(3,3,3));
% 
% Franco
%
% (c) 2012 Stanford VISTA Team

% check inputs
if ~isequal(size(vals,2),size(vol,4)) 
    error('The number measuremnts in vals (%i) is different than the last dimension of the volume (%i)',size(vals,2),size(vol,4))
end
if ~isequal(size(vals,1),size(coords,1)) 
    error('The number values (%i) is different than the number of voxels (%i)',size(vals,2),size(vol,4))
end

% assign the values at the coordinates in the volume.
for iv = 1:length(vals)
    x = coords(iv,1);
    y = coords(iv,2);
    z = coords(iv,3);
    vol(x,y,z,:) = vals(iv,:); 
end

return