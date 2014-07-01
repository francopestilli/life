function feOpenLocalCluster
%
% Initializes a local matlab cluster if the parallel matlab toolbox is
% available.
%
%  feOpenLocalCluster
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com.

% Initialize a local matlab cluster to speed up some of the processes.
if exist('matlabpool','file')
if (matlabpool('size') == 0), 
    if exist('parcluster','file')
    c = parcluster;
    c.NumWorkers = 12;
    else
    matlabpool open;     
    end
else
    disp('[feOpenLocalCluster] Found Matlab parallel cluster open, not intializing.')
end
    disp('[feOpenLocalCluster] Cannot find the Matlab parallel toolbox, not intializing a cluster.')
    disp('[feOpenLocalCluster] Many computations will be substantially slower, without the parallel toolbox.')
end