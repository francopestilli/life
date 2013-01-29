%% s_fe_test_connection
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

recompute_fe = 0; % if 1, it will load a fibergroup and recompute the fe structure.

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
fgFileName    = fullfile(saveDir,'fg','right_hemisphere_occipital_MRTRIX_FG.mat');
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

% Identify the indices of the fibers to remove.
% The ones that do not touch V1
fibersToKeep     = find(~all_fibers);

% Reduce the connectome. We remove all the fibers touching the ROIs, whcih will keep the optic radiation.
feNOTor = feConnectomeSelectFibers(fe,fibersToKeep);

% Now cross-validate the quality fo fit and install the result in the fe structure
fefit   = feFitModel(feGet(feNOTor,'Mfiber'),feGet(feNOTor,'dSig demeaned'),'sgdnn');
feNOTor = feSet(feNOTor,'fit',fefit);
% Save out the fiber group
fgWrite(feGet(feNOTor,'fg acpc'),fullfile(saveDir,'fg','NOT_OR_fg'),'mat');


%% Make an ROI with the voxels of the removed Optic Raadiation fibers
connectome        = feGet(fe,'fg img');
fgOR             = connectome;
fgOR.fibers      = {}; % Clear some memory
fgOR.pathwayInfo = {}; % CLear some memory
OR_fibersToKeep      = find(all_fibers);
for ii = 1:length(OR_fibersToKeep)
 fgOR.fibers{ii,1} = connectome.fibers{OR_fibersToKeep(ii),1};  
end
fgORacpc = dtiXformFiberCoords(fgOR,feGet(fe,'xform img 2acpc'));
fgWrite(fgORacpc,fullfile(saveDir,'fg','optic_radiation_fg'),'mat');

% Compute the value along the optic radiation
% fgORacpcSummary = dtiComputeSuperFiberRepresentation(fgORacpc, [], 100,[]);

% Make an ROI and write it to disk
name      = 'all-OR-fibers';
randColor = rand(1,3);
OR_roi    = dtiNewRoi(name,randColor,fefgGet(fgOR,'unique image coords'));
OR_roiAcpc= OR_roi;
OR_roiAcpc.coords = mrAnatXformCoords(feGet(fe,'xform img2acpc'),OR_roiAcpc.coords);
dtiWriteRoi(OR_roiAcpc, fullfile(saveDir,'rois','optic_radiation_roi'),'mat')

%% Get the RMSE in the volume of the optic radiation-fibers before and after removing the fibers:
% With
rmse_with_or    = (feGet(fe,     'vox rmse',unique(OR_roi.coords,'rows')));

% Save it to file:
niftiName = fullfile(feGet(fe,'savedir'),'rmse_with_or');
rmse_w_map = feValues2volume(rmse_with_or,unique(OR_roi.coords,'rows'),feGet(fe,'map size'));
nii         = niftiCreate('data',rmse_w_map,...
                  'qto_xyz',feGet(fe,'xform img 2 acpc'),...
                  'fname',niftiName);
niftiWrite(nii);

% Without
rmse_without_or = (feGet(feNOTor,'vox rmse',unique(OR_roi.coords,'rows')));

% Save it to file
niftiName   = fullfile(feGet(fe,'savedir'),'rmse_without_or');
rmse_wo_map = feValues2volume(rmse_without_or,unique(OR_roi.coords,'rows'),feGet(fe,'map size'));
nii         = niftiCreate('data',rmse_wo_map,...
                  'qto_xyz',feGet(fe,'xform img 2 acpc'),...
                  'fname',niftiName);
niftiWrite(nii);

%% Get the RMSE in the whole connectome volume before and after removing the fibers:
% With
roi          = feGet(fe,'roi');
rmse_with_or = (feGet(fe, 'vox rmse',unique(roi.coords,'rows')));

% data rmse file
rmseFileName = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/raw/signal_rmse_b1000.nii.gz';


% Save it to file:
[sz]  = dwiGet(dwiFile,'volume size');
niftiName  = fullfile(feGet(fe,'savedir'),'rmse_with_or_wholebrain');
rmse_w_map = feValues2volume(rmse_with_or,unique(roi.coords,'rows'),sz(1:3));
nii        = niftiCreate('data',rmse_w_map,...
                         'qto_xyz',feGet(fe,'xform img 2 acpc'),...
                         'fname',niftiName);
niftiWrite(nii);

% Without
rmse_without_or = (feGet(feNOTor, 'vox rmse',unique(roi.coords,'rows')));

% Save it to file
niftiName   = fullfile(feGet(fe,'savedir'),'rmse_without_or_wholebrain');
rmse_wo_map = feValues2volume(rmse_without_or,unique(roi.coords,'rows'),sz(1:3));
nii         = niftiCreate('data',rmse_wo_map,...
                  'qto_xyz',feGet(fe,'xform img 2 acpc'),...
                  'fname',niftiName);
niftiWrite(nii);


%% make a scatter plot of the loss in explained variance.
mrvNewGraphWin('Increase in RMSE without optic radiation fibers')
plot(rmse_with_or,rmse_without_or,'ro','MarkerFaceColor','r')
axlim = get(gca,'xLim');axis square;set(gca,'yLim',axlim);
hold on, 
plot([axlim(1);axlim(2)], [axlim(1); axlim(2)],'k-')
ylabel('RMSE without optic radiation fibers');xlabel('RMSE with optic radiation fibers')

%% make a histogram RMSE
mrvNewGraphWin('RMSE with and without optic radiation fibers')
[y,x] = hist(rmse_without_or,200);bar(x,y,'r','BarWidth',0.65)
hold on
[y,x] = hist(rmse_with_or,200);bar(x,y,'k','BarWidth',0.35)
ylabel('Occurrences');xlabel('RMSE');
legend({'With optic radiation fibers','Without optic radiation fibers'})

keyboard
%% End