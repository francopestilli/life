function fe = v_lifeCulling
%
% Illustrate how to open up data and and reduce the size of a connectome by
% keeping all the fibers the allow to reduce the size of the connectome
% without reducing the percent variance explained.
% 
% v_lifeCulling
%
%
% Franco (C) 2012 Stanford VISTA team.

%% Initialize a connectome
baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
dtFile     = fullfile(baseDir,'dti40','dt6.mat');
dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
fgFileName = fullfile(baseDir,'fibers','leftArcuateSmall.pdb');

%% Initialize the Connectome
fe = feConnectomeInit(dwiFile,dtFile,fgFileName);

%% Now reduce the size of the fiber groups
% Keep all the fibers that allow not to loose the percent variance
% explained.
fe = feConnectomeCull(fe);
return

%% Save it
feConnectomeSave(fe);

%% Save the fiber group
feConnectomeWrite(fe);

keyboard
