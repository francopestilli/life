function ht_writeFiberDensityNifti(fg, t1, skernal, outname)
% This function will read a fiber group and write a nifti image of fiber
% density in each voxel
%
% ht_writeFiberDensityNifti(fg, t1)
%
% Inputs:
% fg      - Path to fiber group
% t1      - Path to t1.nii.gz
% skernal - Smoothing kernal to apply to the fiber endpoint image
% outname - Full path and file name to save output image

if ~exist('skernal','var') || isempty(skernal)
    skernal = [3 3 3];
end
% Read in the fibers
f = dtiLoadFiberGroup(fg);
% Read in the t1
im = readFileNifti(t1);
% Create an image of fiber density where each voxel counts the number of
% fiber endpoints
fd = dtiComputeFiberDensityNoGUI(f,im.qto_xyz,size(im.data),0,1,1);
% Smooth the image with a guassian kernal 
fd = smooth3(fd,'gaussian',skernal);
% Now build a nifti image that is coregistered to the t1
fdImg = im;
fdImg.data = fd;
fdImg.fname = outname;
% Write the nifti image
writeFileNifti(fdImg);

% see the fiber density on cortical mesh from mrVista window
% use Xform: itkGray(nifti) --> Gray to transform fiber density map (nii
% file) into mrMesh