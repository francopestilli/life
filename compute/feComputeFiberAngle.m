function [az, el, azVec, elVec] = feComputeFiberAngle(fibers)
% Calculate the angle of a fiber at each node.
%
%  [az, el, azVec, elVec] = feComputeFiberAngle(fibers)
%
% INPUTS:
%   fibers  - is a cell array of x,y,z coordinates for each node in a fiber.
%
% OUTPUT:

%   az     - is a cell array of the same length as the number of fibers.
%            Each az{ii} contains a matrix of (numNodes x 1) of azimuth in 
%            degrees.  
%   el     - is a cell array of the same length as the number of fibers.
%            Each az{ii} contains a matrix of (numNodes x 1) of elevation
%            in degrees.      
%
%   azVec  - is a vector of the same length as the number of
%            fibersXnumNodes. It contains all the azimuths for all fibers
%            and nodes in the fiber group.
%   elVec  - is a vector of the same length as the number of
%            fibersXnumNodes. It contains all the elevations for all fibers
%            and nodes in the fiber group.     
% 
% EXAMPLE:
%  fgImg.Pho = feComputeFiberAngle(fibers)
%
% See also: feConnectomeInit.m v_lifeExample.m
%
% Written by Franco Pestilli (c) Stanford University 2013

% Preallocate
nFibers = length(fibers);  % The number of Fibers.
xyz     = cell(1,nFibers); % Memory for the cartesian coordinates of vectors.
az =  cell(1,nFibers);el=az; % Memory for azimuth and elevation in degrees

% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end
parfor ii = 1:nFibers
 % Compute the diffusion gradient at each node of the fiber.
 xyz{ii} = gradient(fibers{ii});
 
 % Now compute the transformation from cartesian coordinates to spherical
 % coordinates.
 [az{ii}, el{ii}] = cart2sph(xyz{ii}(1,:),xyz{ii}(2,:),xyz{ii}(3,:));
 
 % Now transform from radians to deegrees:
 az{ii} = fromRadians('degrees',az{ii});
 el{ii} = fromRadians('degrees',el{ii});
end 

% Cast all azimuths into one big array
azVec = horzcat(az{:});

% Cast all elevations into one big array
elVec = horzcat(el{:});

if ~poolwasopen, matlabpool close; end

return

