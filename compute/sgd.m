function [w, R2] = sgd(y,X,numtoselect,finalstepsize,convergencecriterion,checkerror,nonneg,alpha,lambda)
% 
% This fucntion is working properly and can be used. But the HELP of this
% function has not been updated yet. Blame Franco. F.P. 2012
% 
% Least-square stochastic gradient-descend fit.
% 
% sgd(y,X,numtoselect,finalstepsize,convergencecriterion,checkerror,nonneg)
% 
% <y> is p x 1 with the data
%
% <X> is p x q with the regressors
%
% <numtoselect> is the number of data points to randomly select on each
% iteration
%
% <finalstepsize> is like 0.05
%
% <convergencecriterion> is [A B C] where A is in (0,1), B is a positive
% integer, and C is number of percentages.
%   
%   We stop if we see a series of max(B,round(A*[current total-error-check
%   number])) total-error-check iterations that do not improve performance
%   on the estimation set, where improvement must be better by at least 1%
%   of the previously marked R^2.
% 
% <checkerror> is the number of iterations between total-error-checks
%
% <nonneg> if set to 1 costrains the solution to be positive
% 
% For reference see : Kay et al. 2008 (Supplemental material)
% 
% <alpha, lambda> ElasticNet parameters (optional, defaults to 1 and 0).
% The ElasticNet is a reguarization and variable selection algorithm. 
% The EN penalty is: 
% 
% (y - X*w).^2) + lambda * sum(alpha * w.^2 + (1-alpha) * abs(w)) 
%
% Such that lambda sets the slope of the additional regularization error
% surface and alpha balances between the L1 and L2 constraints. When alpha
% is 1, the algorithm reduces to ridge regression. When alpha is 0, the
% algorithm reduces to the Lasso.
% 
% Reference: Zou and Hastie (Zou & Hastie, (2005) Journal of the Royal
% Statistical Society B 67, Part 2, pp. 301-320)
% See also: Friedman, Hastie and Tibshirani (2008). The elements of
% statistical learning, chapter 3 (page 31, equation 3.54)
%
%
% Copyright Franco Pestilli and Kendrick Kay (2013) Vistasoft Stanford University.

% Set the default for input params: 
if notDefined('convergencecriterion'), convergencecriterion = [.15 3 5]; end 
if notDefined('numtoselect'), numtoselect = 0.1 * size(y,1); end
if notDefined('checkerror'), checkerror=40; end
if notDefined('finalstepsize'), finalstepsize=0.05; end
if notDefined('nonneg'), nonneg = false; end
if notDefined('coordDescent'), coordDescent = false; end

% Set default values for ElasticNet (defaults to regular OLS):
if notDefined('alpha'), alpha   = 0; end 
if notDefined('lambda'), lambda = 0; end

p = size(y,1);  % number of data points
q = size(X,2);  % number of parameters
orig_ssq = full(sum((y).^2)); % Sum of Squres fo the data

% initialize various variables used in the fitting: 
w          = 0 .* rand(q,1); % the set of weights, between 0 and .1
w_best     = w;          % The best set of weights.
est_ssq_best = inf;      % minimum estimation error found so far
estbadcnt  = 0;          % number of times estimation error has gone up
iter       = 1;          % the iteration number
cnt        = 1;          % the total-error-check number
  
% report
fprintf('[%s] Performing fit | %d measurements | %d parameters | ',mfilename,p,q);

% Start computing the fit.
while 1
  % Indices to selected signal and model
  ix      = randi(p,1,numtoselect); % Slower indexing method
  ix2     = false(p,1);
  ix2(ix) = true;

  % select the subset of signal and model to use for fitting
  y0 = y(ix2);
  X0 = X(ix2,:);
  
  % if not the first iteration, adjust parameters
  if iter ~= 1
    % Calculate the gradient (change in error over change in parameter): 
    grad = -((y0 - X0*w)' * X0)' + lambda * (alpha + 2*(1 - alpha)*w);
    
    % This computes the coordinate descent instead of the gradient descent.
    if coordDescent
        % Coordinate descent
        m    = min(grad);
        grad = (grad==m)*m;
    end
    
    % Unit-length normalize the gradient
    grad = unitlengthfast(grad);
    
    % Perform gradient descent
    w = w - finalstepsize*grad;
        
    % Non-negative constrain, we set negative weights to zero
    if ( nonneg ), w(w<0) = 0;end
  end
  
  % check the total error every so often
  if mod(iter,checkerror) == 1
    % Curent estimated sum of the squares of the residuals (SSQ)
    est_ssq = sum((y - X*w).^2) + lambda*(alpha*sum(w) + (1-alpha)*sum(w.^2));
     
    % Check if the SSQ improved
    isimprove = est_ssq < est_ssq_best; 
            
    % We keep fitting if the SSQ is not Inf OR some percent smaller than
    % the best SSQ obtained so far
    %keepfitting = isinf(est_ssq_best) | (est_ssq < ((est_ssq_best - min_ssq)));
    keepfitting = isinf(est_ssq_best) | (est_ssq < ((est_ssq_best * (1-convergencecriterion(3)/100))));

    % do we consider this iteration to be the best yet?
    if isimprove
      % The SSQ was smaller, the fit improved.
      w_best       = w;       % Set the current to be the best so far
      est_ssq_best = est_ssq; % The min error
      
      % OK we improved, but check whether improvement is too small to be
      % considered useful.
      if keepfitting
        % THe fit improved more than the minimum accptable improvement.
        % Reset the counter fo rthe bad fits, so that we start over
        % checking for stopping.
        estbadcnt  = 0;
        %est_ssq_best = est_ssq; % 
      else
        estbadcnt = estbadcnt + 1;
      end
    else
      % The fit actually was bad, SSQ increases count how many bad fit we had. Stop after a centrain number.
      estbadcnt = estbadcnt + 1;
    end
    
    % stop if we haven't improved in a while
    if estbadcnt >= max(convergencecriterion(2),round(convergencecriterion(1)*cnt))
      R2 = 100*(1-(est_ssq_best/orig_ssq));
      fprintf(' DONE fitting | SSQ=%2.3f (Original SSQ=%2.3f) | Rzero-squared %2.3f%%.\n',...
              est_ssq_best,orig_ssq,R2);
      break;
    end
    
    % Update the counter
    cnt = cnt + 1;
  end
  iter = iter + 1;
end

% prepare output
w = w_best;

function [v,len] = unitlengthfast(v,dim)

% function [v,len] = unitlengthfast(v,dim)
%
% <v> is a vector (row or column) or a 2D matrix
% <dim> (optional) is dimension along which vectors are oriented.
%   if not supplied, assume that <v> is a row or column vector.
%
% unit-length normalize <v>.  aside from input flexibility,
% the difference between this function and unitlength.m is that
% we do not deal with NaNs (i.e. we assume <v> does not have NaNs),
% and if a vector has 0 length, it becomes all NaNs.
%
% we also return <len> which is the original vector length of <v>.
% when <dim> is not supplied, <len> is a scalar.  when <dim> is
% supplied, <len> is the same dimensions as <v> except collapsed
% along <dim>.
%
% note some weird cases:
%   unitlengthfast([]) is [].
%   unitlengthfast([0 0]) is [NaN NaN].
%
% example:
% a = [3 0];
% isequalwithequalnans(unitlengthfast(a),[1 0])

if nargin==1
  len = sqrt(v(:).'*v(:));
  v = v / len;
else
  if dim==1
    len = sqrt(sum(v.^2,1));
    v = v ./ repmat(len,[size(v,1) 1]);  % like this for speed.  maybe use the indexing trick to speed up even more??
  else
    len = sqrt(sum(v.^2,2));
    v = v ./ repmat(len,[1 size(v,2)]);
  end
end
