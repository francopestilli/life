function nodesToKeep = feClipFiberNodes(fiber,coords, maxVolDistance)
%
% Clip a fiber to be constrained inside a predefiend volume (x,y,z
% coordinates).
% 
%  function nodesToKeep = feClipFiberNodes(fiber,coords, maxVolDistance)
% 
% INPUTS:
%        fiber          - A fiber, defined as a 3xN array of x,y,z coordinates. 
%                         E.g., fg.fibers{1}.
%        coords         - A volume of defined as a 3xN array of x,y,z coordinates.
%        maxVolDistance - The farther distance (in mm) from the volume a node can
%                         be to be kept in the fiber.
%        
% OUTPUTS:
%        nodesToKeep    - Vector of ones and zeros indicating the nodes to
%                         be kept in the fiber. 1 means "keep the node."
%
% SEE ALSO: feClipFibersToVolume.m, feConnectomePreprocess.m
%
% Franco (c) 2012 Stanford Vista Team.

% Compute the squared distance between each node on fiber ii and the
% nearest roi coordinate
[~, nodesSQdist] = nearpoints(fiber, coords');

% Take the square root of this distance and keep the nodes that are that
% distance inside the volume.
nodesToKeep = sqrt(nodesSQdist) <= maxVolDistance;    

end % Function
