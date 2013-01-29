function [fg xform] = mctNfgStrandCoords2dwi(fg,voxNum)
% 
% function fg = mctNfgStrandCoords2dwi(fg,voxNum)
%
% NFG returns the fibers in an abstract coordinates frame, [-1,1]
%
% here we transform the fibers coodrdinates in the coordinate frame of the
% dwi nifti volume of the sphere they were created from.
%
% voxNum is the number in the volume.
% it can be found in: 8_mri_sim_params.txt
%
% Franco
%
% See also: s_mct_fact_cc_prefrontal_rois.m, mctComputeDataReliability,
%           mctDiffusionModel.m
%
% Example:
%  strands = mctNfgLoadStrands('dir/to/strand/collection');
%  fg = mctNfgStrand2fiber(strands);
%  fg = mctNfgStrandCoords2dwi(fg,voxNum);
% 
% (c) Stanford VISTA Team

if any(cellfun('max',fg.fibers) > 1) || any(cellfun('min',fg.fibers) < -1) 
    error('[%s] the fiber group does not have coordinates between -1 and 1.',mfilename)
end

% this xform trasnforms coordinates from [-1,1] to [0,1]
xform1 = [.5 0 0 .5; 0 .5 0 .5; 0 0 .5 .5];

% this second xform changes coordinates from [0,1] to [0 size(niftiImage)] 
% plus recenters them at in the middle of each voxel
shift = 0.5;    % we shift the coordinates to the center of the voxels
scale = voxNum; % we translate 

% this xform trasforms from [0,1] to image coordinates [0,n] and centers
% the coordinates in the middle of the voxels
xform  = [scale 0 0 shift; 0 scale 0 shift; 0 0 scale shift ];

% now change the coordinates of ach fibers
for ll = 1:numel(fg.fibers)
    % first rescale coordinates between 0 and 1
    fg.fibers{ll} = mrAnatXform(xform1,fg.fibers{ll});
    
    % then scale to dwi image size
    fg.fibers{ll} = mrAnatXform(xform, fg.fibers{ll});
end

return

