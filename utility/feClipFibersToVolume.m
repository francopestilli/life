function fg = feClipFibersToVolume(fg,coords,maxVolDist)
%
% Clip fibers to be constrained withint a Volume.
%
%  function fibers = feClipFibersToVolume(fibers,coords,maxVolDist)
%
% INPUTS:
%        fibers         - A cell array of fibers, each defined as a 3xN array
%                         of x,y,z coordinates. E.g., fg.fibers.
%        coords         - A volume of defined as a 3xN array of x,y,z coordinates.
%        maxVolDist     - The farther distance (in mm) from the volume a node can
%                         be to be kept in the fiber.
%        
% OUTPUTS:
%        fibersOut      - A cell-array of fibers clipped within the volume 
%                        defined by coords.
%
% SEE ALSO: feClipFiberNodes.m, feConnectomePreprocess.m
% 
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

% Handling parallel processing
poolwasopen = 1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end
fibers    = fg.fibers;

% Make sure the coordinates were passed with the expected dimesnions
if ~((size(coords,1) == 3) && (size(coords,2)>=3))
  coords = coords';
end

parfor ii = 1:length(fibers)
  % Compute the squared distance between each node on fiber ii and the
  % nearest roi coordinate
  [~, nodesSqDistance] = nearpoints(fibers{ii}, coords);
  
  % Keep the nodes in the fiber less than the wanted distance
  nodesToKeep = (sqrt(nodesSqDistance) <= maxVolDist);
  fibers{ii}  = fibers{ii}(:,nodesToKeep);
end
fg.fibers     = fibers; clear fibers;

% Remove the empty fibers, fibers that had no nodes in the volume
fibersToKeep = (~cellfun('isempty',fg.fibers));
fg           = fgExtract(fg,fibersToKeep,'keep');

if ~poolwasopen, matlabpool close; end

end % Main function
