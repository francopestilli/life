function pho = feComputeFiberAngle(fibers)
% Calculate the angle of a fiber at each node.
%
%  Q = feComputeFiberAngle(fibers)
%
% INPUTS:
%   fibers  - is a cell array of x,y,z coordinates for each node in a fiber.
%
% OUTPUT:
%   pho     - is a cell array of the same length as the number of fibers.
%             Each pho{ii} contains a matrix of (numNodes x 2) angles, 
%             elevation and azimuth.
%
%
% EXAMPLE:
%  fgImg.Pho = feComputeFiberAngle(fibers,dParms)
%
% See also: feConnectomeInit.m v_lifeExample.m
%
% Franco (c) 2013 Stanford VISTA Team
keyboard % WORK UNFINISHED

% Preallocate
nFibers = length(fibers); % The number of Fibers.
pho       = cell(1,nFibers); % Memory for the angles of each fiber.

% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end
parfor ii = 1:nFibers
 % Compute the diffusion gradient at each node of the fiber.
 fiberGradient = gradient(fibers{ii});
 pho{ii} = fiberGradient;
end 

% Plot on example fiber
fiberIndex = 4;
for fi = 1:length(fibers{fiberIndex})                                                                    
quiver3(fibers{fiberIndex}(1,fi),fibers{fiberIndex}(2,fi),fibers{fiberIndex}(3,fi), ...
    pho{fiberIndex}(1,fi),pho{fiberIndex}(2,fi),pho{fiberIndex}(3,fi));
hold on
end

% Plot the vectors in a signle voxel
% u = sort(unique(fibers(:)))


if ~poolwasopen, matlabpool close; end

return

