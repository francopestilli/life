function fName = feConnectomeSave(fe,varargin)
% Save the 'fe' strucutre.  
%
%   fName = feConnectomeSave(fe,varargin)
%
% INPUTS:
%  savedir - Directory where to save the connectome. 
%            Defualt feGet(fe,'savedir')
%  name    - Name of the connectome.
%            Default feGet(fe,'name')
%
% OUTPUT: fName - path to the save file.
%
% See also, feConnectomeBuild.m, feConnectomeInit.m
%
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

if ~isempty(feGet(fe,'savedir')) && ~exist(feGet(fe,'savedir'),'dir'), 
  mkdir(feGet(fe,'savedir')); 
end

% Get the full file name
fName = feGet(fe,'name');
if ~isempty(varargin) % Add some info to the file name.
    fName = fullfile(feGet(fe,'savedir'),[fName,varargin{1}]);
else
    fName = fullfile(feGet(fe,'savedir'),fName);
end
fprintf('[%s], Saving connectome: %s...\n',mfilename,fName)

% Add a chck here that if a file ith the smae name exists in the directory we append a time tag.
% Otherwise appendign a time tag is not necessary and makes these files less clean.
if (exist([fName,'.mat'],'file') ||  exist(fName,'file'))
  save(fName,'fe','-v7.3','-append');
else
  save(fName,'fe','-v7.3');
end

return