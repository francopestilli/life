function fg_file = feConnectomeWrite(fe,saveDir,fileType)
%
% Writes the Connectome (the fiber group, fe.fg) taking care of the coordinate system. 
%
% fg_file = feConnectomeWrite(fe,[saveDir],[fileType])
%
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

% Default directory for saving the fiber group is the LIFE folder where the
% the fe structure is saved by default.
if notDefined('saveDir'),saveDir = feGet(fe,'save dir');end
if notDefined('fileType'),fileType = 'mat';end

% Build a file name.
fg_file = fullfile(saveDir,sprintf('%s_FG',feGet(fe,'name')));

% Write the fiber group always in acpc coordinates, 
% Also, we always save as .mat
fprintf('[%s] Saving connectome: %s.%s\n',mfilename,fg_file,fileType)
fgWrite(feGet(fe,'fg acpc'),fg_file,fileType); 
% Note .pdb seem to not save correctly by fgWrite. 
% Maybe a problem with pathwaysInfo?

return