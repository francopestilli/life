function feSaveMapToNifti(fe,mapType, niftiName)
%
% Saves a parameter map to file from an fe structure.
%
% feSaveMapToNifti(fe,mapType, niftiName)
%
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

% Generate the requested map and save it to volume.
map = feValues2volume((feGetRep(fe, mapType)),feGet(fe,'roi coords'),feGetRep(fe,'map size'));

% Generate a nifti structure for it.
nii         = niftiCreate('data',map,...
                  'qto_xyz',feGet(fe,'xform img 2 acpc'),...
                  'fname',niftiName);
                
% Write the nifti
niftiWrite(nii);
end
