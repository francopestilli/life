function fe = s_fe_cullingExample(hemisphere)
%
% Illustrate how to open up data and and reduce the size of a connectome by
% keeping all the fibers the allow to reduce the size of the connectome
% without reducing the percent variance explained.
% 
% s_fe_cullinExample
%
%
% Franco (C) 2012 Stanford VISTA team.

dataRootPath = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
subfolders   = fullfile('150dirs_b1000_1');
baseDir      = fullfile(dataRootPath,subfolders);
dtFile       = fullfile(baseDir,'dt6.mat');
dwiFile      = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
switch hemisphere
  case {'right','left'}
    fgFileName   = fullfile(baseDir,'fibers','whole_brain_MRTRIX',sprintf('%s_hemisphere_MRTRIX.mat',hemisphere));
  case {'both','whole'}
    fgFileName   = fullfile(baseDir,'fibers','whole_brain_MRTRIX','WholeBrainFG_MRTRIX.mat');
  otherwise
    fgFileName = hemisphere; % a path to a file name was passed in
end

%% Initialize the Connectome
[~,f,~] = fileparts(fgFileName);
if strcmpi(f(end-1:end),'fe')
   % Load the fe structure
  fprintf('[%s] Loading the fe structure from file: %s...\n',mfilename,fgFileName);
  load(fgFileName);
else % We do this only if we did not get as input an fe structure already
  fe = feConnectomeInit(dwiFile,dtFile,fgFileName);
end

%% Now reduce the size of the fiber groups
% Keep all the fibers that allow not to loose the percent variance
% explained.
fe = feConnectomeCull(fe);

keyboard
return

%% Save it
feConnectomeSave(fe);

%% Save the fiber group
feConnectomeWrite(fe);
