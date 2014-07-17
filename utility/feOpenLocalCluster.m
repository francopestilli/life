function feOpenLocalCluster
% Initializes a local matlab cluster if the parallel matlab toolbox is
% available.
%
%  feOpenLocalCluster
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com.

% Initialize a local matlab cluster to speed up some of the processes.
if exist('matlabpool','file')
   try 
     if (matlabpool('size') == 0)
        if (exist('parcluster','file') == 2)
           c = parcluster;
           c.NumWorkers = 12;
           t = tempname;
           OK = mkdir(t);
           if OK
              c.JobStorageLocation = t;
           end
           matlabpool(c);
        else
           matlabpool open;     
        end
     else
        % disp('[feOpenLocalCluster] Found Matlab parallel cluster open, not intializing.')
     end
   catch ME
     fprintf('\n[feOpenLocalCluster] Problem intializing the cluster: \n\n %s.\n', ME.message)
     disp('[feOpenLocalCluster] Many computations will be substantially slower, without the parallel toolbox.')
   end
else
      disp('[feOpenLocalCluster] Cannot find the Matlab parallel toolbox, not intializing a cluster.')
      disp('[feOpenLocalCluster] Many computations will be substantially slower, without the parallel toolbox.')
end

end
