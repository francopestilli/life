function [fe fefit fexval] = feConnectomeCull(fe,maxNumInter, lowerPrct, percentReduxRMSE)
%
% Find the fibers in a connectome that explain most of the variance, by
% iteratively finding the fibers with largest weights and fitting the life
% model only to those.
%
%   [fe fefit fexval] = feCullConnectome(fe,maxNumInter, lowerPrct, percentReduxRMSE)
% 
% Inputs:
%   fe            - An fe structure, see feCreate.m or v_lifeExample.m
%   maxNumInter - The number of times the fit and fiber selection
%                   processes are repeated.
%   lowerPrct     - The lower percentile for accepting a weight. On each
%                   iteration only fibers above this percentile are
%                   accepted.
%   percentReduxRMSE   - The max percent in R2 reduction that we allow.
%                   When the cross validated R2 get smaller than tis
%                   percent of the original R2, we stop culling.
%
% Outputs: 
%   fe            - The fe structure with  the culled fibers.
%
% Example:
%   baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
%   dtFile     = fullfile(baseDir,'dti40','dt6.mat');
%   dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
%   fgFileName = fullfile(baseDir,'fibers','leftArcuateSmall.pdb');
%   fe         = feConnectomeInit(dwiFile,dtFile,fgFileName);
%   fe         = feConnectomeCull(fe);
% 
% Franco (c) 2012 Stanford VISTA Team.

% For large Conectomes this number might need to be increased, convergence
% might require more interations.
if notDefined('maxNumInter'), maxNumInter = 1000;end

% This is the percent of fibers that are removed at each iteration. The
% larger this number the faster the convergence, but the more likely to
% delete important fibers and produce a reduced connectome with a larger
% loss of R2.
if notDefined('lowerPrct'), lowerPrct     = 10;end

% The percent change in R2 from that of the orignal model. 
% The smaller this number the faster the convergence, the more the fibers
% kept in the connectome.
if notDefined('percentReduxRMSE'), percentReduxRMSE = 2;end

% Fit the model and then reduced it by only accepting fibers that pass
% the minWeights threshold.
for iter = 1:maxNumInter
  fefit  = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');  
   
  % Check whether we start loosing percent variance explained when
  % crossvalidating
  fexval = feXvalidate(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
  
  % Start keeping trak of the quality of fit.
  % When the fit quality decreases we will stop culling.
  if (iter == 1)
    originalRMSE = fexval.rmse.mean;  
    currentRMSE  = fexval.rmse.mean;
  else
    currentRMSE  = fexval.rmse.mean;
  end
  
  if currentRMSE > ((originalRMSE) + (originalRMSE*(percentReduxRMSE/100)))
    break
  end
  fprintf('\n\n[%s] n iter: %i, Original RMSE: %2.3f, Current RMSE: %2.3f.\n\n',mfilename, iter,originalRMSE,currentRMSE)
  
  % Compute the cut-off for the fibers
  minWeight    = prctile(fefit.weights,lowerPrct);
  fibersToKeep = find(fefit.weights > minWeight);
  
  % Reduce the connectome.
  fe = feConnectomeSelectFibers(fe,fibersToKeep);
  
end

% Refit the model and save the fit 
fefit  = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');  
fexval = feXvalidate(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');

% Remove any fit or cross-validation
fe = feSet(fe,'fit', fefit);
fe = feSet(fe,'xvalfit',fexval);

fprintf('[%s] Done culling, in %i iterations.\n Original RMSE: %2.3f, Current RMSE: %2.3f.\n',mfilename, iter,originalRMSE,currentRMSE)

return