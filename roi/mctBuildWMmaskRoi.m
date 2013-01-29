function [roi roiLeft roiRight] = mctBuildWMmaskRoi(dt,faThresh)
% 
% function [roi roiLeft roiRight] = mctBuildWMmaskRoi(dt6,faThresh)
%
% Builds an ROI with FA values above a certain faThreshold.
% 
% The ROI can be used as a white-matter mask before doing other
% computations on the data (e.g. tractography).
%
% Example: 
%    faThresh = 0.01;
%    dt       = dtiLoadDt6('/path/to/dt6/file/');
%    wmMask   = mctBuildWMmask(dt,faThresh); 
%
% See also, test_mictrotrack_simulated_signal_recover
%
% (C) 2012 Stanford VISTA team. 

% if no fa threshold is passed in, 
% return all voxels
if notDefined('faThresh'),faThresh = 0;end

% Compute the FA and normalize
fa = dtiComputeFA(dt.dt6);

% clip bad FA values.
fa(fa>1) = 1; 
fa(fa<0) = 0;

% find the voxels that pass the FA threshold
mask 	 = fa >= faThresh;
[x,y,z]  = ind2sub(size(mask), find(mask));

% Create a new whole-brain ROI from the mask
roi        = dtiNewRoi(sprintf('mask, min fa%2.2f',faThresh));
roi.coords = mrAnatXformCoords(dt.xformToAcpc, [x,y,z]);% Check xfrom if it is correct.

% Left hemisphere, return an roi for all voxels in the left hemisphere
roiLeft = dtiRoiClip(roi, [1 80]);

% Left hemisphere, return an roi for all voxels in the right hemisphere
roiRight = dtiRoiClip(roi, [-80 -1]);

return