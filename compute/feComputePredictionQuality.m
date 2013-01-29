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
% Franco (c) 2012 Stanford VISTA Team

if notDefined('whichr2'), whichr2 = 'percent variance explained'; end

whichr2 = mrvParamFormat(whichr2);

% (1) Compute the quality across voxels.
% R2: coefficient of determination
switch whichr2
  case 'correlation'
    % corrcoef returns a 2x2 and the (1,2) and (2,1) terms are the same and
    % what we want. The (1,1) and (2,2) terms are both 1.
    r2 = 100*(corrcoef(predSig,dSig).^2);
    r2 = r2(1,2);
  case {'coefficientofdetermination','cod'}
    % From Kendrick, should be checked against default percent variance
    % explained, below
    r2 = calccod(dSig,predSig);  % Coefficient of determination
  case {'rzerosquared'}
    % Percent variance explained relative to the data
    % R2        = 1 - (SSE / sum(dSig^2))
    % where SSE = sum[(dSig - predSig).^2]
    r2 = 100 * ...
      (1 - ...
      (sum( (dSig - predSig).^2 ) / ...
      sum( (repmat(mean(dSig),size(dSig,1),1)).^2 ) ) ...
      );
    
  otherwise
    % Percent variance explained
    % R2        = 1 - (SSE / VAR(demeaned dSig))
    % where SSE = sum[(dSig - predSig).^2]
    r2 = 100 * ...
      (1 - ...
      (sum( (dSig - predSig).^2 ) / ...
      sum( (dSig - repmat(mean(dSig),size(dSig,1),1)).^2 ) ) ...
      );
end

% rmse - root mean squared error, this is the most valuable measure,
% although not normalized. it is valuable because the microtrack fits were
% done as a least-squared solution, which assumes that the distribution of
% the residuals and of the rmse is not far from being normally distributed.
rmse = sqrt(mean((dSig - predSig).^2));

return
