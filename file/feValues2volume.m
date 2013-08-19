function vol3d = feValues2volume(vals,coords,vSize)
% Create a 3D volume of sz and place vals at the coord positions
%
%  vol3d = feValues2volume(stat,coords,vSize)
%
% save the vals in the correct locations inside the image volume.
%
% vals   -  N values to insert 
% coords -  N x 3 matrix of coords, (not in acpc)
% vSize  -  The size of the image volume 
% 
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

% Check that we have a value for each coordinate.
if ~( length(vals) == size(coords,1) )
  error('vals (%d) and N coords mis-match (%d)', length(vals),size(coords,1))
end

% Check that the coordinates are contained inside the volume
if all( vSize(1:3) < max(coords) ) || all(min(coords) < [1 1 1])
  error('Roi coordinates are out of range')
end

% Initialize the volume to nan's
vol3d = nan(vSize); 

if (length(vSize) == 3) % This is a statistics that needs to be saved in 3D volume
  for iv = 1:length(vals)
    vol3d(coords(iv,1),coords(iv,2),coords(iv,3)) = vals(iv);
  end
  
elseif (length(vSize) == 4) % This is signal that needs to be saved in 4D volume
  for iv = 1:length(vals)
    vol3d(coords(iv,1),coords(iv,2),coords(iv,3),:) = vals(:,iv)';
  end
  
else
  error('[%s] Volume dimensions mismatch',mfilename)
end