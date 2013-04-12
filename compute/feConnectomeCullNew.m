function [fe, o, fefitGood] = feConnectomeCullNew(fe,maxNumInter, fitType, redux)
%
% Find the fibers in a connectome that explain most of the variance, by
% iteratively finding the fibers with largest weights and fitting the life
% model only to those.
%
%   [fe, o, fefit] = feCullConnectome(fe,maxNumInter, lowerPrct, percentReduxRMSE)
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
% Copyright Franco Pestilli (c) 2013, VISTA Team, Stanford University.

% For large Conectomes this number might need to be increased, convergence
% might require more interations.
if notDefined('maxNumInter'), maxNumInter = 100;end

% The percent of fibers to remove at each iteration.
if notDefined('redux')
    redux.percentRmseIncrease = -2;     % The max percent increase in RMSE from the initial
    redux.percentile = [32 16 8 1]; % The percentile reduction we perform 
    redux.proportionFibers = [0.4, 0.5, 0.6, 0.7]; % The proportion of fiber reduction in which we allow fro each percentile
end

% This is the minimum weight that we use to "keep" fibers.
% We want to delete fibers that have some contribution to the signal.
minWeight    = 0;

% This is the RMSE of the conenctome model as it was passed in. When
% culling fibers, hereafter, we do not want to increase above a certain
% percent of this value.
o.rmseOriginal  = median(feGetRep(fe,'vox rmse'));
o.rmseThreshold = o.rmseOriginal + (redux.percentRmseIncrease*o.rmseOriginal/100);

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
if ~isempty(fe.life.fit)
    fefit = fe.life.fit;   
else
    fefit = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
end
o.weights{1}   = feGet(fe,'fiber weights');
o.numFibers(1) = feGet(fe,'nfibers');
stopped = 0;

% Fit the model and then reduced it by only accepting fibers that pass
% the minWeights threshold.
for iter = 1:maxNumInter
    if iter > 1
        % Re fit the model after eliminating the zero-weighted fascicles
        fefit     = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
        fe        = feSet(fe,'fit', fefit);
    end
    
    % Get the cross-validated RMSE
    o.r2(iter)      = median(feGet(fe,'vox r2'));
    o.r2xv(iter)    = median(feGetRep(fe,'vox r2'));
    o.rmse(iter)    = median(feGet(fe,'vox rmse'));
    o.rmsexv(iter)  = median(feGetRep(fe,'vox rmse'));
    o.rrmse(iter)   = median(feGetRep(fe,'vox rmse ratio'));
    o.weights{iter} = feGet(fe,'fiber weights');
    o.numFibers(iter)    = length(o.weights{iter});
    fprintf('[%s] n iter: %i, Original RMSE: %2.3f, Current RMSE: %2.3f (RMSE Thr %2.3f).\n', ...
        mfilename, iter,o.rmseOriginal, o.rmse(iter),o.rmseThreshold)
    
    % Check whether in the last fit we hit the threshold and increased the
    % rmse. If we did we want to stop, we will keep the fit from the
    % previous iteration not the oen just performed. See below.
    if (o.rmse(iter) > o.rmseThreshold) && iter > 1
        fprintf('[%s] Exiting becuase the RMSE is increasing from initial one...\n',mfilename)
        stopped = 1;
        break  
    end

    %--------------------------------%
    % Now we compute the fibers to remove in two ways.
    % (1) We select a certain percentile of fibers to remove.
    % (2) We select all the fibers with zero weights to remove.
    % (3) We use the method that returns the largest number of removed
    %     fibers.
    
    % (1) Compute percentile:
    %     Sort the weights and get the idices from the sorting
    [~,indx]  = sort(o.weights{iter});
    
    % Find the N percentile and use it to index into the array. This drops
    % only the N percentile.
    % The percent reduction that we apply depends on how many fibers we
    % have already removed. At the beginning of the process we reduce many
    % fibers, as we go along the process we reduce less and less.
    fascicleReduction = 1-o.numFibers(iter)/o.numFibers(1);
    if fascicleReduction < redux.proportionFibers(1)
        whichRedux = 1;
    elseif fascicleReduction < redux.proportionFibers(2)
        whichRedux = 2;
    elseif fascicleReduction < redux.proportionFibers(3)
        whichRedux = 3;
    elseif fascicleReduction < redux.proportionFibers(4)
        whichRedux = 4;        
    end
    fprintf('[%s] PC REDUX fibers to remove: %2.2f%%, actual %2.2f%%, percentile %2.2f\n', ...
            mfilename,redux.proportionFibers(whichRedux), fascicleReduction, ...
            redux.percentile(whichRedux))
    percentileIndex = floor(o.numFibers(iter)*redux.percentile(whichRedux)/100);
    fibersToKeepPC  = indx(percentileIndex:end);
    ftk = false(size(indx));
    ftk(fibersToKeepPC) = true;
    fibersToKeepPC      = ftk; 
    clear ftk

    % (2) Compute the number of zero-weight fibers
    %  Compute the cut-off for the fibers
    fibersToKeepW  = (o.weights{iter} > minWeight);

    % (3) Select the method that returns the largest percentile
    if sum(fibersToKeepPC) < sum(fibersToKeepW)
        fibersToKeep = fibersToKeepPC;
        fprintf('[%s] Using the percentile fiber reduction method (keepPC %i/keepW %i)\n', ...
            mfilename,sum(fibersToKeepPC),sum(fibersToKeepW))
    else
        fibersToKeep = fibersToKeepW;
        fprintf('[%s] Using the 0-weigths fiber reduction method (keepPC %i/keepW %i)\n', ...
            mfilename,sum(fibersToKeepPC),sum(fibersToKeepW))
    end
    %--------------------------------%
    
    % Store the number of fibers removed
    o.removeFibers(iter) = length(find(~fibersToKeep));
    
    % Select the indices fo the fibers that were deleted in the previous
    % loop. The way we address these indices depends on the type of operation.
    o.fibersKept{iter+1} = o.fibersKept{iter}(fibersToKeep);
    if iter > 1
        o.results(iter).r = fefit.results;
    end
    
    fprintf('[%s] n iter: %i, %i current fibers, deleting %i fibers.\n', ...
        mfilename, iter,o.numFibers(iter),o.removeFibers(iter))
   
    if iter > 1
       feGood = fe;
    end
    
    % Reduce the connectome. Meaning remove the fibers that we do not want.
    fe = feConnectomeSelectFibers(fe,find(fibersToKeep));
end

% Save the indices of the fibers that survived all the operations.
o.fibersToKeep( o.fibersKept{end} ) = true;

if stopped
    % Re fit the model after eliminating the zero-weighted fascicles
    fe = feGood;
    fefit     = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
    fe        = feSet(fe,'fit', fefit);
end

% Make sure the size of the M matrix and the weights match, this indicates
% that we have done the due diligence and the fit in the connectome and the
% weights are in sync.
assert(size(feGet(fe,'Mfiber'),2)==size(feGet(fe,'fiber weights'),1));

fprintf('[%s] Done culling, in %i iterations. Original RMSE: %2.3f, Current RMSE: %2.3f.\n',mfilename, iter,o.rmseOriginal,o.rmse(iter))

return
