function test_sgd

%
% Test the stochastic gradient descent algorithm
%
%
% (c) 2012 Stanford VISTA Team


% Set up the regression:
beta1 = rand(10,1);  % All positive!
X1 = randn(1000,10);  % Not necessarily positive:

y = X1 * beta1;


% We should be able to get back the right answer for this simple case:
beta_hat1 = stochasticGD(y,X1, [], [], [], [], true); % With non-negativity constraint
assertAlmostEqual(beta1, beta_hat1, 0.1); % Pretty lenient 

beta2 = randn(10,1);  % Not necessarily positive
X2 = randn(1000,10);  % Not necessarily positive

y = X2 * beta2;

% We should be able to get back the right answer for this simple case:
beta_hat2 = stochasticGD(y,X2);  % No non-negativity constraint this time
assertAlmostEqual(beta2, beta_hat2, 0.1); 

% Sparse one: 
beta3 = [rand(5,1), zeros(5,1)]; 
beta_hat3 = stochasticGD(y,X2,[],[],[],[],[],0.5,1);  % ElasticNet
assertAlmostEqual(beta2, beta_hat2, 0.1); 
