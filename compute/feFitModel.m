function [fit w R2] = feFitModel(M,dSig,fitMethod)
% Fit the LiFE model.
%
% Finds the weights for each fiber to best predict the directional
% diffusion signal (dSig)
%
%  fit = mctDiffusionModelFit(M,dSig,fitMethod)
%
% dSig:  The diffusion weighted signa measured at each
%        voxel in each direction. These are extracted from 
%        the dwi data at some roi coordinates.
% M:     The microtrack difusion model matrix, constructed
%        by feConnectomeBuildModel.m
%
% fitMethod: it can be set to 'lsqnonneg' or 'stochastic gradient descent',
% 'sgd'
%
% See also: feCreate.m, feConnectomeBuildModel.m, feGet.m, feGet.
%
% Example:
%  See v_lifeExample.m
%
% Franco (c) 2012 Stanford VISTA Team

% ** Notes **
%
% The rows of the M matrix are nVoxels*nBvecs. We are going to predict the
% diffusion signal in each voxel for each direction.
%
% The columns of the M matrix are nFibers + nVoxels.  The diffusion signal
% for each voxel is predicted as the weighted sum of predictions from each
% fibers that passes through a voxel plus an isotropic (CSF) term.
%
% In addition to M, we typically return dSig, which is the signal measured
% at each voxel in each direction.  These are extracted from the dwi data
% and knowledge of the roiCoords.
%

% fit the model, by selecting the proper toolbox.
switch fitMethod
  case {'lsqnonneg'}
    options      = optimset('lsqnonneg');
    %options.TolX = '5*eps*norm(c,1)*length(c)';
    w = lsqnonneg(M,dSig,options);
    
  case {'sgd','sgdnn'}% stochastic gradient descend, or non-negative stochastic gradient descend
    tic
    % Stochastic gradient descent method.
    % it solves an L2 minimization problem with non-negative constrain.
    %
    % Basically it takes 'chuncks' of rows of the M matrix and solves those
    % separately but contraining to obtain a consistent global solution.
    signalSiz = size(M,1);
    if signalSiz >= 1000000
      siz     = floor(signalSiz * .1); % size of the chuncks (number rows) taken at every iteration of the solver
    elseif signalSiz > 5000 || signalSiz < 1000000
      siz     = floor(signalSiz * .2); % size of the chunks (number rows) taken at every iteration of the solver
    elseif signalSiz <= 5000
      siz     = signalSiz; % size of the chuncks (number rows) taken at every iteration of the solver
    else
      keyboard
    end
    stepSiz      = 0.04; % step in the direction of the gradient, the larger the more prone to local minima
    stopCriteria = [.15 3 5]; % Stop signals:
    % First, if total error has not decreased less than
    %        an XXX proportion of XXXX.
    % Second, number of small partial fits before
    %         evaluating the quality of the large fit.
    % Third, Amount of R2 improvement judged to be
    %        useful.
    %        It used to be:  percent improvement in R2
    %        that is considered a change in quality
    %        of fit, e.g., 1=1%.
    n      = 50;       % Number of iteration after which to check for total error.
    nonneg = strcmpi(fitMethod(end-2:end),'dnn');
    fprintf('\nLiFE: Computing least-square minimization with Stochastic Gradient Descent...\n')
    [w R2] = sgd(dSig,M,siz, stepSiz, stopCriteria, n,nonneg);
    
    % Save out the Stochastic Gradient Descent parameters
    fit.params.stepSiz      = stepSiz;
    fit.params.stopCriteria = stopCriteria;
    fit.params.numInters    = n;
    
    % Save the state of the random generator so that the stochasit cfit can be recomputed.
    defaultStream = RandStream.getDefaultStream;
    fit.randState = defaultStream.State;   
    
    % Save out some results 
    fit.results.R2        = R2;
    fit.results.nParams   = size(M,2);
    fit.results.nMeasures = size(M,1);
    fprintf(' ...fit process completed in %2.3fs\n',toc)

  otherwise
    error('Cannot fit LiFE model using method: %s.\n',fitMethod);
end

% Save output structure.
fit.weights             = w;
fit.params.fitMethod    = fitMethod;

end
