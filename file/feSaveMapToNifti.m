function feSaveMapToNifti(fe,mapType, niftiName)
% Saves a parameter map to file from an fe structure.
%
% feSaveMapToNifti(fe,mapType, niftiName)
%
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com.

% Generate the requested map and save it to volume.
map = feValues2volume((feGetRep(fe, mapType)),feGet(fe,'roi coords'),feGetRep(fe,'map size'));

% Generate a nifti structure for it.
nii         = niftiCreate('data',map,...
                  'qto_xyz',feGet(fe,'xform img 2 acpc'),...
                  'fname',niftiName);
                
% Write the nifti
niftiWrite(nii);
end
