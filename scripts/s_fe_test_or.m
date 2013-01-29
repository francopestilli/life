%% s_fe_test_or
%
% Illustrate how to compare the quality of predictions between a connectome
% that has a full complement of fascicles, and one that has a
% specific subset of these fascicles removed.
%
% In this example, the full connectome is the occipital lobe.     (C1) We
% then remove all the fibers connected to peri-calcarine (optic radiation). (C2) We use
% LiFE to predict the diffusion signal using both C1 and C2. 
% 
% In a separate script, we visualize the results.  The basic result is that
% in the peri-calcarine region, removing the fibers connected to optic radiation is a
% bad idea. The RMSE of the prediction increases.  It increases mainly in
% that region.
%
% In subsequent scripts we both visualize and create graphs summarizing the
% increased prediction error within the key region of interest.
%
% See also:
%
% Franco (C) 2012 Stanford VISTA team.

%% Start loading the data
%
% At some point, we will put these fe structures into VISTADATA.
% For now, you can compute them, leave them in the feGet(fe,'savedir'), and load them
% up for testing.

recompute_fe = 1; % if 1, it will load a fibergroup and recompute the fe structure.

% Basic directories
if ispc,  basePath = fullfile('\\red.stanford.edu','biac2-wandell6');
else      basePath = fullfile('/biac2','wandell6');
end
dataRootPath  = fullfile(basePath,'data','frk','life_dti','FP20120420');
subfolders    = fullfile('150dirs_b1000_1');
baseDir       = fullfile(dataRootPath,subfolders);
saveDir       = fullfile(baseDir,'LiFE');
dtFile        = fullfile(baseDir,'dt6.mat');
dwiFile       = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
fgFileName    = fullfile(saveDir,'WholeBrainFG_MRTRIX_preprocessed.mat');
roisFileNames = {fullfile(saveDir,'rois','optic_radiation_AND_roi.mat'), ...
                 fullfile(saveDir,'rois','optic_radiation_NOT_roi.mat'), ...
                 fullfile(saveDir,'rois','RV1_useme.mat')};
roiOperations = {'and','not','and'}; % Operations to apply with each ROI

%% Load the data part
if recompute_fe
  % Recompute the fe structure from the original fiber group
  % Initialize the Connectome
  fe = feConnectomeInit(dwiFile,dtFile,fgFileName);
  
  % Fit the model
  fefit  = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsig demeaned'),'sgdnn');
  fe     = feSet(fe,'fit',fefit);
  feConnectomeSave(fe)
  feConnectomeWrite(fe)
  
else
  % Load the fe structure from a precomputed file.
  disp('Loading a precomputed fe structure...')
  feFileName = 'test_culling_more_right_hemisphere_occipital_MRTRIX_FE.mat';
  %'20120830T160637-right_hemisphere_occipital_MRTRIX_FE.mat';
  load(fullfile(saveDir,feFileName))
end

%% Restrict the fibers to the optic radiation
all_fibers = true(feGet(fe,'nfibers'),1);
for iR = 1:length(roisFileNames)
  % Load the ROI
  roi = dtiReadRoi(roisFileNames{iR});
  
  % xform the ROI coordinates to IMG space
  roi.coords =  mrAnatXformCoords(feGet(fe,'xform acpc2img'),roi.coords);
  
  % Keep fibers touching optic radiation
  [~,~, this_fibers] = dtiIntersectFibersWithRoi([], {roiOperations{iR}}, [], roi, feGet(fe,'fg img'));
 
  % Restrict the fiber indices
  all_fibers = and(all_fibers,this_fibers);
end

%% Generate a connectome with all the fibers NOT touching the Optic
% Radiation
feNOTor = feConnectomeSelectFibers(fe,find(~all_fibers));

%% Generate a connectome with ONLY the optic radiation fibers
feOR = feConnectomeSelectFibers(fe,find(all_fibers));

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

% Save memory
clear fe;

%% Use the Optic Radiation connectome as signal for the rest of the model.
% i.e., find the weights of the fibers that expalin the optic radiation in
% the model that should contain NO Optic Radiation
w_OR   = feGet(feNOTor,'Mfiber')\feGet(feOR,'Mfiber');

%% Show the weights of the fit
mrvNewGraphWin('Fit of the mean OR signal by the NOT model');
[y,x] = hist(w_OR(:),ceil(length(w_OR)/50));
bar(x,y,.65);
title(sprintf('Weights\nfor each fiber left in the connectome\nafter removing the optic radiation'));
ylabel('Number of occurrences'), xlabel('Weight value')

%% Plot some drescriptive stats for the weights
mrvNewGraphWin('Weights info');
nTotW  = length(w_OR(:));
minW   = 0.01; % The minimum weight we accept from a fiber, higher-weight fibers are deleted
goodW  = w_OR(:) < minW;
nGoodW = length(find( goodW(:)));
nBadW  = length(find(~goodW(:)));
bar([1 2 3],[nTotW nGoodW nBadW],.65,'k');
set(gca,'box','off','xticklabel',...
  {sprintf('Total (%i)', nTotW),...
   sprintf('Good (%i - %i%%)', nGoodW, floor(100*nGoodW/nTotW)),...
   sprintf('Bad (%i - %i%%)', nBadW,   ceil(100*nBadW/nTotW))},'yscale','log')
title(sprintf('Weights (min w %0.2f)',minW));
ylabel('Number of occurrences'), xlabel('Weight type')


%% Compute the correlation between OR fibers and the Model
ORnorm     = unitlengthfast(feGet(feOR,'Mfiber'));
tic,ORnorm = feGet(feOR,'Mfiber') ./ repmat(sqrt(sum(feGet(feOR,'Mfiber').^2,1)),[size(feGet(feOR,'Mfiber'),1) 1]);  toc % like this for speed.  maybe use the indexing trick to speed up even more??

NOTnorm = unitlengthfast(feGet(feNOTor,'Mfiber'));

%% Now predict the OR fibers with the fitted model
symORfibers = feGet(feNOTor,'M fiber')*w_OR;
prctError   = 100*(symORfibers - feGet(feOR,'M fiber'))./symORfibers;

%% Plot the percent error in fitting 
mrvNewGraphWin('Fit of the OR fibers with the NOT-OR model');
plot(prctError,'-')
title(sprintf('Percent error in fitting.\nBig number means that the removed OR fibers could not be\nexplained by the rest of the fibers.'));
ylabel('Error'), xlabel('Voxel/directions')


keyboard
%% End