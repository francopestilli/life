function fe = v_lifeConnectomeModelProperties
%
% Shwo how to characterize the autocorrelation properties of a Connectome.
% 
% v_lifeConnectomeModelProperties
%
% Franco (C) 2012 Stanford VISTA team.


%% Initialize a connectome
baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
dtFile     = fullfile(baseDir,'dti40','dt6.mat');
dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
fgFileName = fullfile(baseDir,'fibers','leftArcuateSmall.pdb');
  
%% Initialize the Connectome
fe     = feConnectomeInit(dwiFile,dtFile,fgFileName);

%% Estimate the weights and install them in the fe structure
fefit  = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
fe     = feSet(fe,'fit',fefit);

%% Singular Value-Decomposition, find the largest number of singular values
%S = svds((feGet(fe,'model')),600); % THis is a slow full-standard way to
%compute it.
%
% For matrices that have a large number of rows, it is slow to compute the
% SVD in the standard way.
%
% But, we can compute the SVD of the (M'*M), these singular values are
% the S^2 (the square of the svd of M)
Ssquare = svds(((feGet(fe,'model'))'*(feGet(fe,'model'))),feGet(fe,'nfibers'));

%% Plot svd
mrvNewGraphWin('Principal component anaylsis of connectome');
plot(100*cumsum(Ssquare)/sum(Ssquare),'ko-'); % THere is no suaring here inside the sum, because S is square already
hold on
plot([0 feGet(fe,'nfibers')],[98 98],'r-')
ylabel('Percent variance explained ')
xlabel('Number of principal components (svd) in M')

%% Check how correlated fiebrs are in the connectome
% We rmove a random proportion of fibers and use them to generate signal
% for the connectome.
%
% Select a subset of fibers to romove from the connectome.
% Identify the indices of the fibers to remove.
fibersToRemove.percent = 10;
fibersToRemove.num     = round(feGet(fe,'nfibers') .* (fibersToRemove.percent/100));
fibersToRemove.indices = sort(randsample(feGet(fe,'nfibers'),fibersToRemove.num));

% Generate a connectome with all the fibers REMOVED
feTest = feConnectomeSelectFibers(fe,fibersToRemove.indices);

% Generate a connectome with ONLY the optic radiation fibers
allTheOtherFibers = ones(feGet(fe,'nfibers') ,1);
allTheOtherFibers(fibersToRemove.indices) = 0;
feRest = feConnectomeSelectFibers(fe,find(allTheOtherFibers));

%% Fit the Rest to predict the Test
w = feGet(feRest,'M fiber')\feGet(feTest,'M fiber');

%% Predict back the Test fibers
testFibers = feGet(feRest,'M fiber')*w;

%% Plot a histogram of the weights
mrvNewGraphWin('Fit of the Test signal with the Rest model');
[y,x] = hist(w',56); bar(x,y,.65);
title(sprintf('Weights\nfor each fiber left in the connectome\nafter removing the test fibers'));
ylabel('Number of occurrences'), xlabel('Weight value')

%% Plot the percent error in fitting 
mrvNewGraphWin('Fit of the Test signal with the Rest model');
plot(100*(testFibers - feGet(feTest,'M fiber'))./testFibers,'-')
title(sprintf('Percent error in fitting.\nBig number means that the removed fibers could not be\nexplained by the rest of the fibers.'));
ylabel('Error'), xlabel('Voxel/directions')

%% Load the occiptal lobe
% disp('Loading a precomputed fe structure...')
% basePath = fullfile('/biac2','wandell6');
% dataRootPath  = fullfile(basePath,'data','frk','life_dti','FP20120420');
% subfolders    = fullfile('150dirs_b1000_1');
% baseDir       = fullfile(dataRootPath,subfolders);
% saveDir       = fullfile(baseDir,'LiFE');
% feFileName = 'test_culling_more_right_hemisphere_occipital_MRTRIX_FE.mat';
% %'20120830T160637-right_hemisphere_occipital_MRTRIX_FE.mat';
% load(fullfile(saveDir,feFileName))


return

%%
% Render a nice image of the fibers using build3DFrame and
% buildSurfaceFromFrame
fePlot(fe,'connectome');



