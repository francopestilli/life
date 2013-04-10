function [p, EmpiricalDiff] = feComputeConnectivity(w,wo)
%
% function [p, EmpiricalDiff] = feComputeConnectivity(w,wo)
%
% Computes the probability of a connection (fascicle) given the measured
% diffusion data.
%
% INPTUS:
%   w  - the rmse error in each voxel in the volume of the fascicle with
%        the fascicle and all the rest of the fascicles passig trhoguht 
%        the same volume.
%   wo - the rmse in each voxel the volume of the fascicle with out the
%        fascile but with all the rest of the fascicles going through 
%        the same volume.
%
%  OUTPUTS:
%   p  - The probability that the connection represented by the fascicle
%        exists given the data.
%   EmpiricalDiff - The difference in median rmse with and with out the
%                   fascicle. This is the difference observed from the 
%                   data and tested against the null hypothesis of no 
%                   difference with and without the fascicle.
%
% Written by Franco Pestilli (c) 2013 Stanford University

% Compute a test of the diference in rmse
% (1) Get the differece in rmse observed empiriclly
EmpiricalDiff = median(wo) - median(w);

% (2) Compute the Null distribution by:
% (2.1) Combine all the rmse from both WITH and WITHOUT.
% (2.2) Compute 10,000 distributions of rmse one for WITH one for WITHOUT
%       with the same numerosity of the initial samples
% (2.3) Compute the difference between the medians of these 10,000
%       distributions.
NullSet     = [w wo];
sizeWith    = length(w);
sizeWithout = length(wo);

nboots = 1000;
nullDistribution = nan(nboots,1);
parfor ibt = 1:nboots
    % Y = RANDSAMPLE(POPULATION,K)
    BootWith    = median(randsample(NullSet,sizeWith));
    BootWithout = median(randsample(NullSet,sizeWithout));   
    nullDistribution(ibt) = BootWithout - BootWith;
end

% (3) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
if EmpiricalDiff > 0
    if max(nullDistribution)<EmpiricalDiff
        p = 100*1/nboots;
    else
        p = sum(nullDistribution(sort(nullDistribution)>EmpiricalDiff));
    end
else
    if min(nullDistribution)>EmpiricalDiff
        p = -(100*1/nboots);
    else
        p = - sum(nullDistribution(sort(nullDistribution) < EmpiricalDiff));
    end
end
%fprintf('[%s] Probability connection by chance < %2.2f%%\n',mfilename,p)

end

