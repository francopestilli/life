function [fe, o, fefitGood] = feConnectomeCullExample(fe,maxNumInter, fitType, percentReduxRMSE)
%
% Find the fibers in a connectome that explain most of the variance, by
% iteratively finding the fibers with largest weights and fitting the life
% model only to those.
%
% NOTE: This function is a derivation of feConnectomeCull.m with the
% intention to show how by reducing too many fibers the quality of fit is
% reduced. It does this by removing a 15% of the fibers at every iteration.
%
%   [fe, o, fefit] = feCullConnectomeExample(fe,maxNumInter, lowerPrct, percentReduxRMSE)
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
% Franco (c) 2013 Stanford VISTA Team.

% For large Conectomes this number might need to be increased, convergence
% might require more interations.
if notDefined('maxNumInter'), maxNumInter = 1000;end

% The percent change in R2 from that of the orignal model. 
% The smaller this number the faster the convergence, the more the fibers
% kept in the connectome.
if notDefined('percentReduxRMSE'), percentReduxRMSE = 10;end

% This is the RMSE of the conenctome model as it was passed in. When
% culling fibers, hereafter, we do not want to increase above a certain
% percent of this value.
o.rmseOriginal  = median(feGetRep(fe,'vox rmse'));
o.rmseThreshold = o.rmseOriginal + (percentReduxRMSE*o.rmseOriginal/100);

% Initialize some outputs:
o.removeFibers = nan(maxNumInter,1);
o.rmse         = nan(maxNumInter,1);
o.rmsexv       = nan(maxNumInter,1);
o.r2           = nan(maxNumInter,1);
o.r2xv         = nan(maxNumInter,1);
o.rrmse        = nan(maxNumInter,1);
o.numFibers    = nan(maxNumInter,1);

% The following is the vector containing the indices 
% to the fibers we KEEP from the origianl fiber group
% after all the logical operations are applied.
o.fibersToKeep = false(feGet(fe,'nfibers'),1); 

% The following is a cell array which will old the relative indices of the
% fibers into each sized-down version of the fibers.
% Each entry of the cell arry holds the indices to the fibergroup in the
% before the current operation was applied.
o.fibersKept{1} = 1:length(o.fibersToKeep);

fefit = fe.life.fit;
lambda = length(feGet(fe,'dsigdemeaned'))*2;
stopped = 0;
% Fit the model and then reduced it by only accepting fibers that pass
% the minWeights threshold.
for iter = 1:maxNumInter
    if iter > 1
        % Re fit the model after eliminating the zero-weighted fascicles
        fefitGood = fefit;
        fefit     = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),fitType,lambda);
        fe        = feSet(fe,'fit', fefit);
    end
    
    % Get the cross-validated RMSE
    o.r2(iter)      = median(feGet(fe,'vox r2'));
    o.r2xv(iter)    = median(feGetRep(fe,'vox r2'));
    o.rmse(iter)    = median(feGet(fe,'vox rmse'));
    o.rmsexv(iter)  = median(feGetRep(fe,'vox rmse'));
    o.rrmse(iter)   = median(feGetRep(fe,'vox rmse ratio'));
    fprintf('[%s] n iter: %i, Original RMSE: %2.3f, Current RMSE: %2.3f (RMSE Thr %2.3f).\n',mfilename, iter,o.rmseOriginal, o.rmse(iter),o.rmseThreshold)
    
    % Compute the cut-off for the fibers
    o.weights{iter}      = feGet(fe,'fiber weights');
    
    % Sort the weights and get the idices from the sorting
    [~,indx]  = sort(o.weights{iter});
    
    % Find the N percentile and use it to index into the array. THis drops
    % only the N percentile.
    percentileIndex = floor(size(o.weights{iter},1)*percentReduxRMSE/100);
    fibersToKeep = indx(percentileIndex:end);

    % Store the number of fibers removed
    o.removeFibers(iter) = percentileIndex;
    o.numFibers(iter)    = length(o.weights{iter});
    
    % Select the indices fo the fibers that were deleted in the previous
    % loop. The way we address these indices depends on the type of operation.
    o.fibersKept{iter+1} = o.fibersKept{iter}(fibersToKeep);
    if iter > 1
        o.results(iter).r = fefit.results;
    end
    
    fprintf('[%s] n iter: %i, %i current fibers, deleting %i fibers.\n',mfilename, iter,o.numFibers(iter),o.removeFibers(iter))
    
    % Reduce the connectome.
    fe = feConnectomeSelectFibers(fe,(fibersToKeep));
end

% Save the indices of the fibers that survived all the operations.
o.fibersToKeep( o.fibersKept{end} ) = true;

fefit     = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn',lambda);
fe        = feSet(fe,'fit', fefit);

% Make sure the size of the M matrix and the weights match
assert(size(feGet(fe,'Mfiber'),2)==size(feGet(fe,'fiber weights'),1));

fprintf('[%s] Done culling, in %i iterations. Original RMSE: %2.3f, Current RMSE: %2.3f.\n',mfilename, iter,o.rmseOriginal,o.rmse(iter))

return
