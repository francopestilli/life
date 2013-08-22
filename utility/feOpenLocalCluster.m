function feOpenLocalCluster
%
% Initializes a local matlab cluster if the parallel matlab toolbox is
% available.
%
%  feOpenLocalCluster
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

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
    disp('[feOpenLocalCluster] Cannot find the Matlab prallel toolbox, nto intializing a cluster.')
    disp('[feOpenLocalCluster] Many computations will be substantially slower, withut the parallel toolbox.')
end

end