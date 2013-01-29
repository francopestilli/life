function ni = feWriteValues2nifti(stats,fileName,xForm2Acpc)
% Write the vals in the nifti file named fileName 
% 
% ni = feWriteValues2nifti(stats,fileName,xForm2Acpc)
%
%
% Write statistics computed for each voxel in an ROI into a nifti image.
% 
% stats is a 3d volume containing values at the coordinate sof an roi.
%
% 
% See also: 
%
% Example:
%    feWriteValues2nifti(fd,'fiberDensity_of_some_roi',xform)
%
% Franco (c) 2012 Stanford VISTA Team 

% Create the nifti structure
ni  = niftiCreate('data',stats, ...
    'qto_xyz',xForm2Acpc, ...
    'fname',fileName, ...
    'nifti_type',class(stats));

% Write the file.  Saved as a nii.gz
niftiWrite(ni)

end

