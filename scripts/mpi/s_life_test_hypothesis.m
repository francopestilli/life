function s_life_test_hypothesis(fitType,save)
% function s_life_test_hypothesis
%
% (1) Loads a fascicle produced with ConTrack, i.e., tracts between two ROIs. 
% (2) Generates a White-Matter ROI as the volume identifed by all the
%     voxels of fascicle.
% (3) Runs life on the ConTrack fibers to select the non-zero weighted
%     fibers.
% (4) Uses mrtrix to track inside the ROI.
% (5) Runs life on the connectome returned by MRTRIX to select the non-zero
%     weighted fibers.
% (6) Runs LiFE and computes R2 again to establish the quality of fit.
% (7) Adds the ConTrack fascicle to the connectome.
% (8) Fits LiFE and computes R2.
% (9) Computes the R2 in (6) with that in (8), which one is better?
%
% Implements a split-in-half cross validation on the same data-set, by
% fitting the model on one half and testing the prediction on the other
% half.
%
% Example: 
%   s_life_test_hypothesis
%
% Franco (C) 2012 Stanford VISTA team.

if notDefined('dtiDataType'); dtiDataType = 'b1000'; end
if notDefined('lambda');           lambda = 0;       end
if notDefined('fitType');         fitType = 'sgd';   end
if notDefined('save');               save = 0;       end % Save or not files and figures
if notDefined('algo');algo = 'MRTRIX';end % STT, TEND or MRTRIX

if ispc
     baseDir  = '\\white.stanford.edu\biac2-wandell6\data\arokem\ModelFits\FP20120420\';
else
    baseDir  = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
end

dtDir    = fullfile('150dirs_b1000_1');
dtFile   = fullfile(baseDir,dtDir,'dt6.mat');
dwiFile  = fullfile(baseDir,'raw',sprintf('0009_01_DWI_2mm150dir_2x_%s_aligned_trilin.nii.gz',dtiDataType));

fiberDir = fullfile(baseDir,'LiFE','optic_radiation_test');
saveDir  = fiberDir;
saveName = sprintf('%s_%s_%s',mfilename,dtiDataType,fitType);

% Set default parameters
[fe dwi] = mctSetDefaultParms(algo,dtiDataType,fiberDir,dtFile,dwiFile,baseDir,fitType,lambda,saveDir,saveName,save);
keyboard
% (1) load the contrack fascicle for the optic radiation
fe = loadFiberGroupLocal(fe);

fprintf('\n[%s] Analyzing Fascicle %s\n',mfilename,fe.fg.img.name)
% (2) Generate a white-matter ROI out of it.
fe.fg.img.coords = fefgGet(fe.fg.img,'unique image coords');

% (3) Run LiFE on the ConTrack fascicle.
fe = mctRunLifeFitsLocal(fe,dwi);

% (4) Track using MRTRIX in the volume identified by the contrack fascicle.

% (5) Run LiFE on the connectome returned by MRTRIX, remove the
%     zero-weighted fibers.
fe = mctRunLifeFitsLocal(fe,dwi);

% (6) Now cross-validate the quality fo fit
fe = mctCrossValidate(fe);

% (7) Add the Contrack fascicle to the MRTRIX connectome.

% (8) Run LiFe to the combined connectome.
fe = mctRunLifeFitsLocal(fe,dwi);

% (9) Now cross-validate the quality fo fit
fe = mctCrossValidate(fe);

% Plot results and fiber groups
fe = mctPlotResultsLocal(fe);

% Save the R2 and rmse per voxel in a nifti file.
mctSaveQualityNifti(fe,dwi)

% Now save a new nifti image with the predicted signal
mctSavePredictionNifti(fe,dwi)

end


%%%%%%%%%%%%%%%%%%%%%%%
% mctSaveQualityNifti %
%%%%%%%%%%%%%%%%%%%%%%%
function mctSaveQualityNifti(fe,dwi)
%
% Takes the predicted signal of a life model and saves it into a new nifti
% file identical to the original DWI file except for the signal where the
% the fiber group made predictions.
%
% Franco

% Saving a dwi nifti image with the predicted signal.
fprintf('[%s] Saving a dwi nifti images with the predicted and residual signals.\n',mfilename)

% zero-out the data field of the nifti first
dwi.nifti.data   = zeros(size(dwi.nifti.data,1),size(dwi.nifti.data,2),size(dwi.nifti.data,3));

% change the nifti size information
dwi.nifti.dim    = size(dwi.nifti.data);
dwi.nifti.pixdim = dwi.nifti.pixdim(1:3);
dwi.nifti.ndim   = 3;

% Get the fascicle name and dwi-filename
fascicleName                          = fe.fg.img.name;
fascicleName(isspace(fe.fg.img.name)) = '_';
[~, dwiName] = fileparts(dwi.nifti.fname);
[~, dwiName] = fileparts(dwiName);

% Save the rmse and R2 of LiFE
% we replace only the directions, not the B0, saved in the first
dwi.nifti.data = mctReplaceImageValues(dwi.nifti.data,fe.life.vox.r2,fe.coords);

% Save the file
fileName = fullfile(fe.savedir, sprintf('LiFE_R2_%s_%s.nii.gz',fascicleName,dwiName));
niftiWrite(dwi.nifti,fileName)

% zero-out the data field of the nifti first
dwi.nifti.data = zeros(size(dwi.nifti.data,1),size(dwi.nifti.data,2),size(dwi.nifti.data,3));
dwi.nifti.data = mctReplaceImageValues(dwi.nifti.data,fe.life.vox.rmse,fe.coords);

% Save the file
fileName = fullfile(fe.savedir, sprintf('LIFE_RMSE_%s_%s.nii.gz',fascicleName,dwiName));
niftiWrite(dwi.nifti,fileName)

% zero-out the data field of the nifti first
dwi.nifti.data = zeros(size(dwi.nifti.data,1),size(dwi.nifti.data,2),size(dwi.nifti.data,3));
dwi.nifti.data = mctReplaceImageValues(dwi.nifti.data,fe.orig.vox.r2,fe.coords);

% Save the file
fileName = fullfile(fe.savedir, sprintf('ORIG_R2_%s_%s.nii.gz',fascicleName,dwiName));
niftiWrite(dwi.nifti,fileName)

% zero-out the data field of the nifti first
dwi.nifti.data = zeros(size(dwi.nifti.data,1),size(dwi.nifti.data,2),size(dwi.nifti.data,3));
dwi.nifti.data = mctReplaceImageValues(dwi.nifti.data,fe.orig.vox.rmse,fe.coords);

% Save the file
fileName = fullfile(fe.savedir, sprintf('ORIG_RMSE_%s_%s.nii.gz',fascicleName,dwiName));
niftiWrite(dwi.nifti,fileName)

end


%%%%%%%%%%%%%%%%%%%%
% mctCrossValidate %
%%%%%%%%%%%%%%%%%%%%
function fe = mctCrossValidate(fe)
%
% We use split in Half cross validation.
%
% This can be improved but calculations become slow if we make the cross
% validation scheme more interesting.
%

% Corss validate the LiFE model
evenHalf = ones(length(fe.dSig_demeaned),1);
evenHalf(2:2:end) = 0;
evenHalf = logical(evenHalf);
oddHalf  = ~evenHalf;

% Fit the LiFE to odd and even halfs
evenW = mctFitDiffusionModel(fe.Afiber(evenHalf,:),fe.dSig_demeaned(evenHalf),'sgd');
oddW  = mctFitDiffusionModel(fe.Afiber(oddHalf,:), fe.dSig_demeaned(oddHalf), 'sgd');

% Predicted signal with the weights generated by the other half of the signal.
evenPSig  = mctComputePredictedSignal(fe.Afiber(evenHalf,:),oddW);
oddPSig   = mctComputePredictedSignal(fe.Afiber(oddHalf,:),evenW);

% Compute quality of fit.
[evenRmse evenR2] = mctComputePredictionQuality(fe.dSig_demeaned(evenHalf),evenPSig);
[oddRmse oddR2]   = mctComputePredictionQuality(fe.dSig_demeaned(oddHalf),oddPSig);

% Output the mean and STD of the cross-validate quality of fit.
fe.life.xval.rmse(1) = mean([evenRmse oddRmse]);
fe.life.xval.rmse(2) = std([evenRmse  oddRmse]);
fe.life.xval.r2(1)   = mean([evenR2   oddR2]);
fe.life.xval.r2(2)   = std([evenR2    oddR2]);

% Cross validate the Original model
Aorig  = sum(fe.Afiber,2); % rebuild the model

% Predicted signal with the weights generated by the other half of the signal.
evenW = Aorig(evenHalf) \ fe.dSig_demeaned(evenHalf);
oddW  = Aorig(oddHalf) \ fe.dSig_demeaned(oddHalf);

evenPSig = mctComputePredictedSignal(Aorig(evenHalf),oddW);
oddPSig  = mctComputePredictedSignal(Aorig(oddHalf),evenW);

[evenRmse evenR2] = mctComputePredictionQuality(fe.dSig_demeaned(evenHalf),evenPSig);
[oddRmse oddR2]   = mctComputePredictionQuality(fe.dSig_demeaned(oddHalf),oddPSig);

% Output the mean and STD of the cross-validate quality of fit.
fe.orig.xval.rmse(1) = mean([evenRmse oddRmse]);
fe.orig.xval.rmse(2) = std([evenRmse  oddRmse]);
fe.orig.xval.r2(1)   = mean([evenR2   oddR2]);
fe.orig.xval.r2(2)   = std([evenR2    oddR2]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%
% mctReplaceImageValues %
%%%%%%%%%%%%%%%%%%%%%%%%%
function img = mctReplaceImageValues(img,vals,coords,indexes)
% 
% img = mctReplaceImageValues(img,vals,coords,indexes)
%
% Takes a 3/4D image (img, like a dwi data file) and 
% replaces the values of the image at the coordinates (coords)
% with the values in vals.
%
% Franco

% Indexes specify the entries in the image where to copy the new values
% for example if the first 10 entries are the B0-values and we want to 
if notDefined('indexes')
  indexes = 1:size(img,4);
end
if ~( length(indexes) <= size(img,4) )
  error('Image and vals do not have the same size.')
end

if ~( size(vals,2) == size(coords,1) )
  error('Vals and coords do not have the same size.')
end

if ~( length(indexes) == size(vals,1) )
  error('Image and vals do not have the same size.')
end

for ic = 1:size(coords,1)
  img(coords(ic,1),coords(ic,2),coords(ic,3),indexes)   = vals(:,ic);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% mctSavePredictionNifti %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function mctSavePredictionNifti(fe,dwi)
%
% Takes the predicted signal of a life model and saves it into a new nifti
% file identical to the original DWI file except for the signal where the
% the fiber group made predictions.
%
% Franco

% Saving a dwi nifti image with the predicted signal.
fprintf('[%s] Saving a dwi nifti images with the predicted and residual signals.\n',mfilename)

% Get and set some basic info.
bvalsIndx = dwiGet(dwi,'b0imagenums');
nImages   = dwiGet(dwi,'n images');
theseImages = (bvalsIndx(end)+1):nImages;

% Get the fascicle name and dwi-filename
fascicleName                          = fe.fg.img.name;
fascicleName(isspace(fe.fg.img.name)) = '_';
[~, dwiName] = fileparts(dwi.nifti.fname);
[~, dwiName] = fileparts(dwiName);

% Save the signal predicted by LiFE
% Reshape the signal to the volume
reshaped_pSig = mctReshape(fe.life.pSigFull, fe.nBvecs, fe.nVoxels);

% Now replace the image values at the ROI coordinates.
% we replace only the directions, not the B0, saved in the first
dwi.nifti.data = mctReplaceImageValues(dwi.nifti.data,reshaped_pSig,fe.coords,theseImages);

% Save the file
fileName = fullfile(fe.savedir, sprintf('LiFE_predicted_%s_%s.nii.gz',fascicleName,dwiName));
niftiWrite(dwi.nifti,fileName)

% Save the residuals of the original model.
% Reshape the signal to the volume
reshaped_res = mctReshape(fe.life.resSig, fe.nBvecs, fe.nVoxels);
% Now replace the image values at the ROI coordinates.
% we replace only the directions, not the B0, saved in the first
dwi.nifti.data = mctReplaceImageValues(dwi.nifti.data,reshaped_res,fe.coords,theseImages);

% Save the file
fileName = fullfile(fe.savedir, sprintf('LiFE_residuals_%s_%s.nii.gz',fascicleName,dwiName));
niftiWrite(dwi.nifti,fileName)

% Save the signal predicted by the original model.
% Reshape the signal to the volume
reshaped_pSig = mctReshape(fe.orig.pSigFull, fe.nBvecs, fe.nVoxels);
% Now replace the image values at the ROI coordinates.
% we replace only the directions, not the B0, saved in the first
dwi.nifti.data = mctReplaceImageValues(dwi.nifti.data,reshaped_pSig,fe.coords,theseImages);

% Save the file
fileName = fullfile(fe.savedir, sprintf('ORIG_predicted_%s_%s.nii.gz',fascicleName,dwiName));
niftiWrite(dwi.nifti,fileName)

% Save the residuals of the original model.
% Reshape the signal to the volume
reshaped_res = mctReshape(fe.orig.resSig, fe.nBvecs, fe.nVoxels);
% Now replace the image values at the ROI coordinates.
% we replace only the directions, not the B0, saved in the first
dwi.nifti.data = mctReplaceImageValues(dwi.nifti.data,reshaped_res,fe.coords,theseImages);

% Save the file
fileName = fullfile(fe.savedir, sprintf('ORIG_residuals_%s_%s.nii.gz',fascicleName,dwiName));
niftiWrite(dwi.nifti,fileName)

end


%%%%%%%%%%%%%%%%%%%%%%
% mctSetDefaultParms %
%%%%%%%%%%%%%%%%%%%%%%
function [fe dwi] = mctSetDefaultParms(algo,dtiDataType,fiberDir,dtFile,dwiFile,baseDir,fitType,lambda,saveDir,saveName, save)
%
% Set the dafault parameters for the two algorithms being tested.
%
% Franco

% load the dwi and dti files        
dwi          = dwiLoad(dwiFile);
[dtiF, dtiH] = mrDiffusion('off',dtFile);

% Set parameters for the TEND fiber group
fe.algo           = algo;
fe.dtiDataType    = dtiDataType;
fe.fg.name        = fullfile(fiberDir,'mct_fg_r_optic_radiation_contrack.mat');
fe.xform.img2acpc = dtiGet(dtiH,'img 2 acpc xform');
fe.xform.acpc2img = dtiGet(dtiH,'acpc 2 img xform');
fe.fitType        = fitType;
fe.lambda         = lambda;
fe.savedir        = [saveDir,'_',fe.algo];
if ~isdir(fe.savedir), mkdir(fe.savedir);end
fe.save           = save; % save figures and results
fe.savename       = saveName;

% Close mrDiffusion
close(dtiF), drawnow

% load the T1, this is used to display a slice with the fiber groups
%disp('Loading T1')
%nifti = readFileNifti(fullfile(baseDir,'t1','t1.nii.gz'));

end


%%%%%%%%%%%%%%%%%%%%%%%
% loadFiberGroupLocal %
%%%%%%%%%%%%%%%%%%%%%%%
function fe = loadFiberGroupLocal(fe)
%
% Load the contrack fibers and the MRTRIX, STT or TEND fibers.
%
fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,fe.fg.name)
fe.fg.acpc = dtiLoadFiberGroup(fe.fg.name);
fe.fg.img  = dtiXformFiberCoords(fe.fg.acpc, fe.xform.acpc2img,'img');

end


%%%%%%%%%%%%%%%%%%%%
% runLifeFitsLocal %
%%%%%%%%%%%%%%%%%%%%
function fe = mctRunLifeFitsLocal(fe,dwi)
%
% Performs a run of the life model and evaluates it against the
% original model
%
% Franco

% Only compute the ROI coordinates for the independent ROI case.
fe.coords  = fefgGet(fe.fg.img, 'uniqueimagecoords');

fprintf('\n[%s] the ROI used has size(%i,%i)\n',mfilename,size(fe.coords))

% Build the LiFE model
[Afiber, Aiso, fe.dSig, fe.dSig_demeaned] = mctBuildDiffusionModel(dwi,fe.fg.img,fe.coords,[],'ones');

% Number of fibers
fe.nFiber = size(Afiber,2); 

% Output the A matrix to use later for cross-validation
fe.Afiber = Afiber;

% build the full model
Afull = [Afiber, Aiso];

% Fit the model.
fe.life.w.fiber = mctFitDiffusionModel(Afiber, fe.dSig_demeaned, fe.fitType,fe.lambda);

% Fit the isotropic
fe.life.w.iso  = Aiso \ fe.dSig;
fe.life.w.full = [fe.life.w.fiber; fe.life.w.iso];

% Predict the signal 
fe.life.pSigFull  = mctComputePredictedSignal(Afull,fe.life.w.full);
fe.life.pSig      = mctComputePredictedSignal(Afiber,fe.life.w.fiber);

% Overall quality across voxels
[fe.life.rmse, fe.life.r2] = mctComputePredictionQuality(fe.dSig_demeaned, fe.life.pSig,1);

% Compute residuals, this is done on the demeaned (fiber) signal. 
% This signal was not explained by the model 
fe.life.resSig = (fe.dSig - fe.life.pSigFull) + (fe.dSig - fe.dSig_demeaned);

% Compute the diffusion signal predicted by the original fiber model.
fe.nVoxels = size(Aiso,2); 
fe.nBvecs  = size(Afiber,1) / fe.nVoxels;
AorigFull  = [sum(Afiber,2)  Aiso * fe.life.w.iso]; 
w          = AorigFull \ fe.dSig;

fe.orig.w.fiber  = w(1);
fe.orig.w.iso    = w(2);
fe.orig.w.full   = w;
fe.orig.pSigFull = mctComputePredictedSignal(AorigFull,fe.orig.w.full);
fe.orig.pSig     = mctComputePredictedSignal(AorigFull(:,1),fe.orig.w.fiber);

% Overall quality across voxels
[fe.orig.rmse, fe.orig.r2] = mctComputePredictionQuality(fe.dSig_demeaned, fe.orig.pSig,2);

% Compute residual signal
fe.orig.resSig = (fe.dSig - fe.orig.pSig) + AorigFull(:,2);

% Reshape the signals by voxel
fe.nBvecs        = dwiGet(dwi,'n diffusion images');
fe.nVoxels       = size(fe.coords,1);
fe.life.vox.pSig = mctReshape(full(fe.life.pSig), fe.nBvecs, fe.nVoxels);
fe.life.vox.dSig = mctReshape(full(fe.dSig_demeaned), fe.nBvecs, fe.nVoxels);

[rmse, r2]       = mctComputePredictionQuality(fe.life.vox.dSig, fe.life.vox.pSig,2);
fe.life.vox.r2   = r2;
fe.life.vox.rmse = rmse;
fe.life.vox.res  = mctReshape(full(fe.life.resSig), fe.nBvecs, fe.nVoxels);

fe.orig.vox.pSig   = mctReshape(full(fe.orig.pSig), fe.nBvecs, fe.nVoxels);
fe.orig.vox.pSigUW = mctReshape(full(fe.orig.pSig) * (1/fe.orig.w.fiber), fe.nBvecs, fe.nVoxels);
fe.orig.vox.dSig   = mctReshape(full(fe.dSig_demeaned), fe.nBvecs, fe.nVoxels);
[rmse, r2]         = mctComputePredictionQuality(fe.orig.vox.dSig, fe.orig.vox.pSig,2);
fe.orig.vox.r2     = r2;
fe.orig.vox.rmse   = rmse;
fe.orig.vox.res    = mctReshape(full(fe.orig.resSig), fe.nBvecs, fe.nVoxels);

% Select the fibers explaining most of the variance
fe.life.w.criterion(1) = 0.75; % get the top 75% fibers
fe.life.w.criterion(2) = quantile(fe.life.w.fiber,fe.life.w.criterion(1));
fe.life.hits           = fe.life.w.fiber >= fe.life.w.criterion(2); % the fibers with large weights
fe.life.fa             = ~fe.life.hits;

end


%%%%%%%%%%%%%%%%%%%%%%%
% mctPlotResultsLocal %
%%%%%%%%%%%%%%%%%%%%%%%
function fe = mctPlotResultsLocal(fe)
%
% plots results of a life fit.
%
% Franco

% plotting fiber predictions against the signal  
fe.figs.hFit  = mctDisplayModelFitLocal(fe);
% plot quality of fit.
fe.figs.hQual = mctPlotQualityOfFit(fe.life.r2,fe.orig.r2, fe.life.pSig,fe.orig.pSig, fe.dSig_demeaned);

if fe.save
  fe = feSaveResultsLocal(fe);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mctDisplayModelFitLocal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = mctDisplayModelFitLocal(fe)
% 
% function h = mctDisplayModelFitLocal(wFibers,wIso,fitType)
%
% Display the microtrack fits, still working on this. 
%
% Franco (C) 2012 Stanford VISTA team. 

h = mrvNewGraphWin(sprintf('LiFE Results - Using: %s',fe.fitType));

% Fiber weights
cLiFE = [.4576 .5424 1];cOrig = [1 .5 0];
subplot(3,2,1),
plot(fe.life.w.fiber,'o','MarkerFaceColor',cLiFE,'Color',cLiFE), axis tight
hold on,
plot([0 length(fe.life.w.fiber)],[fe.orig.w.fiber fe.orig.w.fiber],'-','Color',cOrig), axis tight
set(gca,'yLim',[-.1 max(fe.life.w.fiber) + 0.01],'TickDir','out','box','off');
title('LiFE solution - Fibers weights')
xlabel('Fiber number');ylabel('Fiber weight')

% isotropic components' weights, mean signal in a voxel
subplot(3,2,2), plot(fe.life.w.iso,'k.'), axis tight
set(gca,'YLim',[-0.1 nanmax(fe.life.w.iso(:))*1.1],'TickDir','out','box','off');
title('Mean DW signal in each voxel')
xlabel('Voxel number / diffusion direction index');ylabel('DW signal')

% Plot the signals in the voxels
subplot(3,2,3:4) 
plot(fe.dSig_demeaned,'k-'); hold on
plot(fe.life.pSig,'-','LineWidth',1,'MarkerFaceColor',cLiFE,'Color',cLiFE,'MarkerEdgeColor','w');
plot(fe.orig.pSig,'--','LineWidth',1,'MarkerFaceColor',cOrig,'Color',cOrig,'MarkerEdgeColor','w');
axis tight
ylabel('Demeaned DW signal')
title(sprintf('Percent variance explained: Orig:%2.2f, LiFE:%2.2f',fe.orig.r2,fe.life.r2))
set(gca,'TickDir','out','box','off');

% Plot the signal in the voxels given by the origianl model before
% optimization
subplot(3,2,5:6) 
plot(fe.orig.pSig * 1/fe.orig.w.fiber,'--','LineWidth',1,'MarkerFaceColor',cOrig,'Color',cOrig,'MarkerEdgeColor','w');
hold on
plot(fe.dSig_demeaned,'k-');
axis tight
xlabel('Voxels and diffusion directions');
ylabel('Demeaned DW signal')
title('Original fiber model and measured data')
set(gca,'TickDir','out','box','off');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mctSaveLifeResultsLocal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fe = feSaveResultsLocal(fe)
%
% Save the results to disk.
%
% Franco

% Remove blank spaces in files name.
fe.fg.acpc.name(isspace(fe.fg.acpc.name))='_';
fe.file.name  = sprintf('%s_%s_%s',fe.savename,fe.algo,fe.fg.acpc.name);
fe.file.results = fullfile(fe.savedir,fe.file.name);
fprintf('\n\nSaving file: %s\n\n',fe.file.results)
save(fe.file.results,'-v7.3')

disp('Saving results figures ... ')
% eps for the plots
% saveas(fe.figs.hFit ,fullfile(fe.savedir,[fe.file.name,'_Fit_',fe.algo]),'fig')
savefigvista(fe.figs.hFit ,[fe.file.name,'_Fit_',fe.algo],'eps',fe.savedir,1,0);

% saveas(fe.figs.hQual ,fullfile(fe.savedir,[fe.file.name,'_Qual_',fe.algo]),'fig')
savefigvista(fe.figs.hQual ,[fe.file.name,'_Qual_',fe.algo],'eps',fe.savedir,1,0);

if isfield(fe.figs,'hAll')
  
  % fig for the fibers
  disp('Saving fibers figures ... ')
  saveas(fe.figs.hAll(1) ,fullfile( fe.savedir,[fe.file.name,'_FgAll_',fe.algo]),'fig')
  saveas(fe.figs.hGood(1) ,fullfile(fe.savedir,[fe.file.name,'_FgGood_',fe.algo]),'fig')
  saveas(fe.figs.hBad(1) ,fullfile( fe.savedir,[fe.file.name,'_FgBad']),'fig')
end

end

