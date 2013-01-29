load /biac4/wandell/biac2/wandell6/data/frk/LiFE/results/mct_rt12_tfocsl2_b1000_lambda0.mat;

Afiber         % 533100x18788
dSig_demeaned  % 533100x1
fiber_w        % 1 x 18788  (nonnegative)

calccod(Afiber*fiber_w',dSig_demeaned)   % 57.134

% 2 hours

%%%%%

h = stochasticGD( ...
  dSig_demeaned, ...
  Afiber, ...
  10000, ...     % bigger means gradient is more robust but potentially slower
  .05, ...       % step size
  [.1 10 1], ... % stop if total error has not decreased recently
  10);           % how often to evaluate total error
calccod(Afiber*h,dSig_demeaned)

%%%%%

% no early stopping [fit all the way!]
% no bootstrapping [ok]
% no missing data [ok]
% no z-scoring issues [HRM]
% no momentum [let's not worry about this speed issue yet]
% no dc estimation [already eliminated]
% no precomputation [let's not worry about this speed issue yet]
% no sparsity [doesn't make sense giving that we are minimizing sq error fully]
% stochastity of the fitting results seems unavoidable.  (different results each time; can't really get to global minimum)

% references: 
% - kay 08b uses gradient descent.  see the supplementary materials.
% - the current code just does a stochastic version of gradient descent (google "stochastic gradient descent")
