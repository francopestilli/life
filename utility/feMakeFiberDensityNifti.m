function fdNii = feMakeFiberDensityNifti(fe, mapName, smoothingKernel,best)
% 
% This function will extract the fiber group from an fe structure 
% and write a nifti image of fiber density in each voxel.
%
%   fdNii = feMakeFiberDensityNifti(fe)
%
% INPUTS:
%   fe              - An fe structure.
%   smoothingKernel - 3D smoothing kernel to apply to the fiber endpoint
%                     image. If set to 0 no smoothing will be applied.
%   mapName         - Full path and file name to save output image
%   best            - Allows to show all the fibers (1) or only the ones that
%                     have non-zero weight (0). Default is to show all the
%                     fibers (0).
%
% OUTPUT:
%   fdNifti -  nifti file containing the fiber density map
%
% Franco (c) Stanford Vista Team, 2013 

% The 3D smoothing kernel to apply to the nifti image so that the density
% will look more uniform
if ~exist('smoothingKernel','var') || isempty(smoothingKernel),
  smoothingKernel = [3 3 3];
end

% Build a file name if it was not passed in.
if ~exist('mapName','var') || isempty(mapName),
  mapName = [feGet(fe,'name'),'_fiberDensity'];
end

% Extract all the fibers
fg = feGet(fe,'fibers acpc');
 
% Select the ones that have non-zero weights.
% If requested, default = no.
if ~exist('best','var') || isempty(best), best= 0;end
   
if best
  w  = feGet(fe,'fiber weights');
  fg.fibers = fg.fibers(w > 0);
end

% Load the high-res anatomical file and to build a nifti image that is 
% coregistered with the segmentation and surface files.
fdNii       = niftiRead(feGet(fe,'anatomy file'));
fdNii.fname = mapName;

% Create an image of fiber density where each voxel counts the number of
% fiber endpoints
fdNii.data = dtiComputeFiberDensityNoGUI(fg, ...
                                 fdNii.qto_xyz,    ...
                                 size(fdNii.data), ...
                                 0,1,1);
clear fg

if all(smoothingKernel) > 0
  % Smooth the image with a guassian kernel
  fdNii.data = smooth3(fdNii.data,'gaussian',smoothingKernel);
end

% Write the nifti image
niftiWrite(fdNii);

% see the fiber density on cortical mesh from mrVista window
% use Xform: itkGray(nifti) --> Gray to transform fiber density map (nii
% file) into mrMesh

return