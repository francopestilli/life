function [meanBestRD meanBestAD medianBestRD medianBestAD stdBestRD stdBestAD bestRD bestAD] = mctComputeTensorParams(dt,faCutoff)
% function [meanBestRD meanBestAD medianBestRD ...
%           medianBestAD stdBestRD stBestAD bestRD bestAD] = mctComputeTensorParams(dt,[faCutoff])
%
% This function takes a dt structure and computes the mean median and std
% for the Axial and Radial diffusivity for the whole brain volume.
%
% These values will be then used to build the tensors for a fiber group
% created in this brain volume.
% 
% Example:
%         [d_ad d_rd] = mctComputeTensorParams(dt);
%        
%         dParms(1) = d_ad; dParms(2) = d_rd; dParms(3) = d_rd;
%         fg.Q = fgTensors(fg,dParms);
%
% (C) 2011 Stanford VISTA team. 


% compute RD and AD
[fa,~,rd,ad] = dtiComputeFA(dt.dt6);

% set the nan's and unfeasable fa values to 0's
fa(isnan(fa))  = 0;
fa(fa >= 1)    = 0;

% sort the best and worse fa
[fasorted, indexes] = sort(fa(:));

% set a min FA value 
if notDefined('faCutoff')
   % get a cutoff fa value
   % 10% most anisotropic voxels 
   temp     = ceil(length(fasorted) - lenght(fasorted)*0.1);
   faCutoff = fasorted(temp);
end

% select the most anisotropic voxels
inx = find(fa > faCutoff);

% find the RD diffusiovity for the most anisotropic voxels
bestRD       = rd(inx);
meanBestRD   = mean(bestRD(:));
medianBestRD = median(bestRD(:));
stdBestRD    = std(bestRD(:));

% find the RD diffusiovity for the most anisotropic voxels
bestAD       = ad(inx);
meanBestAD   = mean(bestAD(:));
medianBestAD = median(bestAD(:));
stdBestAD    = std(bestAD(:));

return