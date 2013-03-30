function [az, el, azVec, elVec, azByVox, elByVox] = feComputeFiberAngle(fg)
% Calculate the angle of a fiber at each node.
%
%  [az, el, azVec, elVec, azByVox, elByVox] = feComputeFiberAngle(fg)
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
% NOTES:
%   Elevation is measured from the x-y plane. Notice that if elevation = 0, 
%   the point is in the x-y plane. If elevation = pi/2, then the point is 
%   on the positive z-axis.
%
% Written by Franco Pestilli (c) Stanford University 2013

% Preallocate
nFibers      = fefgGet(fg,'nFibers');% The number of Fibers.
xyz     = cell(1,nFibers); % Memory for the cartesian coordinates of vectors.
az =  cell(1,nFibers);el=az; % Memory for azimuth and elevation in degrees
% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

% Get the indices of the nodes in eahc voxel
allVoxels = fefgGet(fg,'uniqueimagecoords');

tmp_fibers = fg.fibers;
parfor ii = 1:nFibers
 % Compute the diffusion gradient at each node of the fiber.
 xyz{ii} = gradient(tmp_fibers{ii});
 
 % Now compute the transformation from cartesian coordinates to spherical
 % coordinates. Azimuth and Elevation.
 [az{ii}, el{ii}] = cart2sph(xyz{ii}(1,:),xyz{ii}(2,:),xyz{ii}(3,:));
 
 % Now transform from radians to degrees: 
 % 
 % Elevation is measured from the x-y plane. Notice that if elevation = 0, 
 % the point is in the x-y plane. If elevation = pi/2, then the point is 
 % on the positive z-axis.
 [az{ii}, el{ii} ] = fromRadians('degrees',az{ii},el{ii});
end 
clear tmp_fibers

% The following is the mapping of fiber nodes to voxels in the roi. 
nodes2voxels = fefgGet(fg,'nodes2voxels',allVoxels);

% Now repack the angles in terms of voxels:
tic
nCoords    = size(allVoxels,1);
voxelsInFG = fefgGet(fg,'voxels in fg',nodes2voxels);
elByVox    = cell(1,nCoords);azByVox=elByVox;

for thisFiber=1:nFibers
    voxelsInFiber = voxelsInFG{thisFiber};   % A few voxels, in a list
    azInFiber  = az{thisFiber}; % The corresponding azimuth in the fiber
    elInFiber  = el{thisFiber}; % The corresponding elevation in the fiber

    % Then add a row for each (fiber,node) pairs that pass through
    % the voxels for this fiber.
    for jj=1:length(voxelsInFiber)
        thisVoxel = voxelsInFiber(jj);
        azByVox{thisVoxel} = cat(1,azByVox{thisVoxel},azInFiber(jj));
        elByVox{thisVoxel} = cat(1,elByVox{thisVoxel},elInFiber(jj));
    end
end
fprintf('[%s] azimuth and elevation in each voxel completed in: %2.3fs.\n',mfilename, toc)
    
% Cast all azimuths into one big array
azVec = horzcat(az{:});

% Cast all elevations into one big array
elVec = horzcat(el{:});

if ~poolwasopen, matlabpool close; end

return

