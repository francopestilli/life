function n = mctComputeRiceNoise(mu,sd,siz)
%
% function n = mctComputeRiceNoise(mu,sd,siz)
%
% Generates Rician distributed noise.
% 
% Rician noise can be added to the diffusion signal to generate simulations.
%
% Example:
%     n = mctComputeRiceNoise(500,10,[10,10]);
%
% See also,
%     s_mct_simulated
%
% Franco
%
% (C) 2012 Stanford VISTA team.

% Matlab does not provide a Rician Noise generator. To approximate rician
% noise I add the square of two gaussian distributions with the same
% standard deviation (sd), one with zero mean the other with some mean
% (mu). The square root operation takes me back to the original space.
%
% See "Related Distribution" here:
%
% http://en.wikipedia.org/wiki/Rice_distribution

% generate using the square root of the sum of the squares of two gaussians
x = mu + sd .* randn(siz);
y =      sd .* randn(siz);
n = sqrt(x.^2 + y.^2);

% The following is an alternative way of approximating Rician noise. In
% general this one generates narrower distributions. But these
% distributions have with longer tails when sd is small compared to mu
% (i.e., when the distribution looks normal, e.g., mu=40, sd=100), they are
% just narrower when sd is large compared to mu (when the distribution does
% not looks normal, e.g., mu=40,sd=10)

%p = random('poiss',mu^2/(2*sd^2));
%x = random('chi2',2*p+2,siz);
%n = sd*sqrt(x);

% Show the generated distribution
%mrvNewGraphWin('Rician Noise');
%hist(n(:),length(n(:))*0.1);
