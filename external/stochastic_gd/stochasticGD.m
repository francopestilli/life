function w = stochasticGD(y,X,numtoselect,finalstepsize,convergencecriterion,checkerror,nonneg,alpha,lambda)
% 
% This fucntion is working properly and can be used. But the HELP of this
% function has not been updated yet. Blame Franco. F.P. 2012
%
%
% Least-square stochastic gradient-descend fit.
% 
% stochasticGD(y,X,numtoselect,finalstepsize,convergencecriterion,checkerror,nonneg)
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
% Kendrick & Franco (c) 2012 Stanford VISTA team.

% Set the default for input params: 
if notDefined('convergencecriterion'), convergencecriterion = [.15 3 5]; end 
if notDefined('numtoselect'), numtoselect = 0.1 * size(y,1); end
if notDefined('checkerror'), checkerror=40; end
if notDefined('finalstepsize'), finalstepsize=0.05; end
if notDefined('nonneg'), nonneg = false; end

% Set default values for ElasticNet (defaults to regular OLS):
if notDefined('alpha'), alpha = 0; end 
if notDefined('lambda'), lambda = 0; end

p     = size(y,1);    % number of data points
q     = size(X,2);    % number of parameters
totsq = sum((y).^2);  % total squared error of the data.

% initialize various variables used in the fitting: 
w          = zeros(q,1); % the set of weights, between 0 and .1
w_best     = w;          % The best set of weights.
%esterr = [];            % record of estimation set error
esterr_min = Inf;        % minimum estimation error found so far
estr2_min  = NaN;        % best R^2 found so far
estbadcnt  = 0;          % number of times estimation error has gone up
iter = 1;                % the iteration number
cnt = 1;                 % the total-error-check number

% report
fprintf('[%s] Performing fit | %d measurements | %d parameters | ',mfilename,p,q);

% do it
while 1
  % indices of data points to select
  %ix = ceil(rand(1,numtoselect)*p);% slower
  ix = randi(p,1,numtoselect);     % faster 2x

  % select
  y0 = y(ix);
  X0 = X(ix,:);
  
  % if not the first iteration, adjust parameters
  if iter ~= 1
    
    % Least-squares solution.
    %
    % calculate gradient (change in error over change in parameter): 
    grad = -((y0 - X0*w)' * X0)' + lambda * (alpha + 2*(1 - alpha)*w);
    
    % unit-length normalize the gradient
    grad = unitlengthfast(grad);
    
    % perform gradient descent
    w = w - finalstepsize*grad;
        
    % Non-negative cosntrain, we set the weights below zero to zero
    if ( nonneg ), w(w<0) = 0;end
  end
  
  % check the total error every so often
  if mod(iter,checkerror) == 1
    
    % Least-squares solution.
    %
    % calculate and record sum of the squares of the residuals
    esterr0 = sum((y - X*w).^2) + lambda*(alpha*sum(w) + (1-alpha)*sum(w.^2));
    
    estr2 = 100*(1-esterr0/totsq);
    %esterr = [esterr esterr0];
    isimprove = esterr0 < esterr_min;  % did we improve sq error
    % % did we actually improve enough over previous checkpoint:
    
    % The covergence criterion here is coded as a percentage of improved
    % R2 from the previous iteration.
    isgood = isnan(estr2_min) | estr2 > estr2_min * (1+convergencecriterion(3)/100);
    
    % The covergence criterion here is coded as the amount of improved R2
    % over the previous iteration. Not the percentage.
    %isgood = isnan(estr2_min) | estr2 > estr2_min + convergencecriterion(3);
    
    % report
    %fprintf('iter: %03d | est: %.3f / R2=%.1f%% (%s)\n',iter,esterr0,estr2,knkChoose(isimprove,'D','U'));
    
    % do we consider this iteration to be the best yet?
    if isimprove
      % numiters = iter;  % the iteration number of the best set of weights so far
      w_best = w;         % the best set of weights so far
      esterr_min = esterr0;
      
      % if we improved enough over previous checkpoint, reset the counter
      if isgood
        estbadcnt = 0;
        estr2_min = estr2;
      else
        estbadcnt = estbadcnt + 1;
      end
    else
      estbadcnt = estbadcnt + 1;
    end
    
    % stop if we haven't improved in a while
    if estbadcnt >= max(convergencecriterion(2),round(convergencecriterion(1)*cnt))
      fprintf(' DONE fitting | R2=%2.3f.\n',(100*(1-esterr_min/totsq)));
      break;
    end
    
    cnt = cnt + 1;
    
  end
  
  iter = iter + 1;
  
end

% prepare output
w = w_best;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUNK:

% VISUALIZE
%     subplot(1,2,1); plot(esterr,'ro-');
%     subplot(1,2,2); bar(w(round(linspace(1,length(w),500))));
%     drawnow;


%     case -1
%       val = thresh*max(abs(grad));
%       grad(grad > -val & grad < val) = 0;
%       grad = unitlengthfast(grad);
%       w = w - finalstepsize*grad;
%       w(w<0) = 0;
%
%     switch fittype
%     case 0
%     case 1
%       [d,ix] = sort(abs(grad),'descend');
%
%       ixcnt = 1;
%       while 1
%         delta = finalstepsize*sign(grad(ix(ixcnt)));
%         if w(ix(ixcnt)) - delta >= 0
%           w(ix(ixcnt)) = w(ix(ixcnt)) - delta;
%           break;
%         else
%           ixcnt = ixcnt + 1;
%         end
%       end
%     end
%
