function  [rmse, r2] = feComputePredictionQuality(dSig,predSig,whichr2)
% 
% function  [rmse, r2] = feComputePredictionQuality(dSig,predSig,whichr2)
%
% Compute the quality of prediction between two data sets.
%
% dste = data standard error, it is a measure of the variability of the measurements
%       or a measure of the quality of the fit. it is a vector of sandard
%       errors of the measurements, the vector has the same length of dSig and predSig. 
%       i.e., dste_j = sqrt(sum(data_i - mean(data_i))^2 / (n - 1)) / sqrt(n)
%       where j = 1:length(dSig) and i = 1:numRepeatedDtiScans 
%
% See also: s_mct_fact_cc_prefrontal_rois.m, mctComputeDataReliability,
%           mctDiffusionModel.m
%
% Example:
%  See [rmse, r2] = feComputePredictionQuality(dSig,predSig)
% 
% (c) Franco Stanford VISTA Team

if notDefined('whichr2'), whichr2 = 1; end

% (1) Compute the quality across voxels.
% R2: coefficient of determination
if whichr2 == 1
  r2 = 100*(corr(predSig,dSig).^2);   % Can this become corrcoef???
elseif whichr2 == 2
  r2 = 100 * (1 - (sum( (dSig - predSig).^2 ) ./ sum( (dSig - repmat(mean(dSig),size(dSig,1),1)).^2 ) ));
else
  r2 = calccod(dSig,predSig);  % Coefficient of determination
end

% rmse - root mean squared error, this is the most valuable measure,
% although not normalized. it is valuable because the microtrack fits were
% done as a least-squared solution, which assumes that the distribution of
% the residuals and of the rmse is not far from being normally distributed.
rmse = sqrt(mean((dSig - predSig).^2));

% (2) Compute quality by voxel.
%dSig    = reshape(dSig,nBvecs,nVoxels);
%predSig = reshape(predSig,nBvecs,nVoxels);

% We compute the the coefficient of determination and root-mean-squared
% error by voxel
%r2vox = coeffDetermination(predSig,dSig);
%rmseVox = sqrt(mean((dSig - predSig).^2));

return
