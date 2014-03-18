function [fit w R2] = feFitModel(M,dSig,fitMethod,lambda)
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
    fprintf('\nLiFE: Computing least-square minimization with LSQNONEG...\n')
    options      = optimset('lsqnonneg');
    w = lsqnonneg(M,dSig,options);
    fprintf(' ...fit process completed in %2.3fs\n',toc)
    R2=[];
  case {'bbnnls'}
    tic
    fprintf('\nLiFE: Computing least-square minimization with BBNNLS...\n')
    opt = solopt;
    opt.maxit = 500;
    opt.use_tolo = 1;
    out_data = bbnnls(M,dSig,zeros(size(M,2),1),opt);
    fprintf('BBNNLS status: %s\nReason: %s\n',out_data.status,out_data.termReason);
    w = out_data.x;
    fprintf(' ...fit process completed in %2.3f hours\n',toc/60/60)
    % Save the state of the random generator so that the stochasit cfit can be recomputed.
    defaultStream = RandStream.getGlobalStream; %RandStream.getDefaultStream;
    fit.randState = defaultStream.State;   
    
    % Save out some results 
    fit.results.R2        = [];
    fit.results.nParams   = size(M,2);
    fit.results.nMeasures = size(M,1);
    R2=[];
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
    elseif signalSiz > 10000 || signalSiz < 1000000
      siz     = floor(signalSiz * .5); % size of the chunks (number rows) taken at every iteration of the solver
    elseif signalSiz <= 10000
      siz     = signalSiz; % size of the chuncks (number rows) taken at every iteration of the solver
    else
      keyboard
    end
    stepSiz      = 0.0124; % step in the direction of the gradient, the larger the more prone to local minima
    stopCriteria = [.1 5 1]; % Stop signals:
    % First, if total error has not decreased less than
    %        an XXX proportion of XXXX.
    % Second, number of small partial fits before
    %         evaluating the quality of the large fit.
    % Third, Amount of R2 improvement judged to be
    %        useful.
    %        It used to be:  percent improvement in R2
    %        that is considered a change in quality
    %        of fit, e.g., 1=1%.
    n      = 100;       % Number of iteration after which to check for total error.
    nonneg = strcmpi(fitMethod(end-2:end),'dnn');
    fprintf('\nLiFE: Computing least-square minimization with Stochastic Gradient Descent...\n')
    [w, R2] = sgd(dSig,M,siz,        stepSiz,      stopCriteria,        n,         nonneg);
             %sgd(y,   X,numtoselect,finalstepsize,convergencecriterion,checkerror,nonneg,alpha,lambda)
    % Save out the Stochastic Gradient Descent parameters
    fit.params.stepSiz      = stepSiz;
    fit.params.stopCriteria = stopCriteria;
    fit.params.numInters    = n;
    
    % Save the state of the random generator so that the stochasit cfit can be recomputed.
    defaultStream = RandStream.getGlobalStream; %RandStream.getDefaultStream;
    fit.randState = defaultStream.State;   
    
    % Save out some results 
    fit.results.R2        = R2;
    fit.results.nParams   = size(M,2);
    fit.results.nMeasures = size(M,1);
    fprintf(' ...fit process completed in %2.3fs\n',toc)

    case {'sgdl1','sgdl1nn'}% stochastic gradient descend, or non-negative stochastic gradient descend
    tic
    % Stochastic gradient descent method.
    % it solves an L2 minimization problem with non-negative constrain.
    %
    % Basically it takes 'chuncks' of rows of the M matrix and solves those
    % separately but contraining to obtain a consistent global solution.
    signalSiz = size(M,1);
    if signalSiz >= 1000000
      siz     = floor(signalSiz * .1); % size of the chunks (number rows) taken at every iteration of the solver
    elseif signalSiz > 10000 || signalSiz < 1000000
      siz     = floor(signalSiz * .5); % size of the chunks (number rows) taken at every iteration of the solver
    elseif signalSiz <= 10000
      siz     = signalSiz; % size of the chunks (number rows) taken at every iteration of the solver
    else
      keyboard
    end
    stepSiz      = 0.0124; % step in the direction of the gradient, the larger the more prone to local minima
    stopCriteria = [.1 5 1]; % Stop signals:
    % First, if total error has not decreased less than
    %        an XXX proportion of XXXX.
    % Second, number of small partial fits before
    %         evaluating the quality of the large fit.
    % Third, Amount of R2 improvement judged to be
    %        useful.
    %        It used to be:  percent improvement in R2
    %        that is considered a change in quality
    %        of fit, e.g., 1=1%.
    n      = 100;       % Number of iteration after which to check for total error.
    nonneg = 1;
    fprintf('\nLiFE: Computing least-square minimization (L1) with Stochastic Gradient Descent...\n')
    %lambda = [length(dSig)*2.75];
    [w, R2] = sgdL1(dSig,M,siz, stepSiz, stopCriteria, n,nonneg,[],lambda);
    fprintf('Lambda: %2.2f | nFibers: %i | L1 penalty: %2.3f | L2 penalty: %2.3f\n',lambda, length(find(w>0)),sum(w),sum(w.^2))

    % Save out the Stochastic Gradient Descent parameters
    fit.params.stepSiz      = stepSiz;
    fit.params.stopCriteria = stopCriteria;
    fit.params.numInters    = n;
    
    % Save the state of the random generator so that the stochasit cfit can be recomputed.
    defaultStream = RandStream.getGlobalStream; %RandStream.getDefaultStream;
    fit.randState = defaultStream.State;   
    
    % Save out some results 
    fit.results.R2        = R2;
    fit.results.nParams   = size(M,2);
    fit.results.nMeasures = size(M,1); 
    fit.results.l2        = sum(w.^2);
    fit.results.l1        = sum(w);
    
    fprintf(' ...fit process completed in %2.3fs\n',toc)

  otherwise
    error('Cannot fit LiFE model using method: %s.\n',fitMethod);
end

% Save output structure.
fit.weights             = w;
fit.params.fitMethod    = fitMethod;

end
