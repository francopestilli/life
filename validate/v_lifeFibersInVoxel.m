function fe = v_lifeFibersInVoxel
%
% Illustrate how to open up data and run a small linear fascicle evaluation
% (LIFE).
% 
% v_lifeFibersInVoxel
%
%
% Written by Franco Pestilli (c) 2013 Stanford VISTA team.


%% Initialize a connectome
baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
dtFile     = fullfile(baseDir,'dti40','dt6.mat');
dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
fgFileName = fullfile(baseDir,'fibers','leftArcuateSmall.pdb');
t1FileName = fullfile(baseDir,'t1','t1.nii.gz');
repFile    = fullfile(baseDir,'raw','dwi.nii.gz');

%% Read the fascicle
fg = fgRead(fgFileName);
allCoords = fefgGet(fg, 'uniqueimagecoords');

% Pick a voxel out of the ROI
vox = 500;
roiCoords = allCoords(vox,:);

%% Clip the connectome to be constrained within the volumeRoi.
tic, fprintf('Clipping Connectome fibers that leave the volume ROI...\n');
fg = feClipFibersToVolume(fg,roiCoords,2);
fprintf('process completed in %2.3fhours\n',toc/60/60);

%% Initialize the Connectome
fe     = feConnectomeInit(dwiFile,dtFile,fg,mfilename,[],repFile,t1FileName);

%% Estimate the weights and install them in the fe structure
fefitL1  = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdl1nn',10000);
fe     = feSet(fe,'fit',fefitL1);

% Now Remove half of the fascicles, the ones with highest weights
[~, keepFascicles] = sort(feGet(fe,'fiber weights'));
keepFascicles = keepFascicles(1:length(keepFascicles)/2);
feRedux       = feConnectomeReduceFibers(fe, keepFascicles );

% Now fit again
feRedux  = feSet(feRedux,'fit',feFitModel(feGet(feRedux,'Mfiber'),feGet(feRedux,'dsigdemeaned'),'sgdl1nn',10000));


return



