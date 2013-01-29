% s_lifeExampleFibersProperties
%
% Illustrate how to open up data, run a small linear fascicle evaluation
% (LIFE) and compute some properties of the fibers.
% 
% s_lifeExampleFibersProperties
%
% Franco (C) 2012 Stanford VISTA team.


%% Initialize a connectome
baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
dtFile     = fullfile(baseDir,'dti40','dt6.mat');
dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
fgFileName = fullfile(baseDir,'fibers','leftArcuateSmall.pdb');

%% Compute a bunch of descriptive poperties of the fibers.
% We compute these for each node in eahc fiber. 
% These values depend on the values in the dtFile. 
fg   = fgRead(fgFileName);
fLen = fefgGet(fg,'length');
eig  = fefgGet(fg,'eigenvals',dtFile);
fa   = fefgGet(fg,'fa',eig);
md   = fefgGet(fg,'md',eig);
rd   = fefgGet(fg,'rd',eig);
ad   = fefgGet(fg,'ad',eig);
wsh  = fefgGet(fg,'westinshape',eig);

return
