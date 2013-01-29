function s_life_to12_test(fitType, lambda,algo, this_roi)
% function s_life_to12_test
%
% Load the TO12 roi, load some tracks produced with STT or TEND out of such
% roi. 
% 
% Then load the same tracks without the ILF.
%
% Then produce a fit to the DW data in the LGN ROI using the original
% tracks summed with the Optic Radiation tracks.
%
% Then produce a fit to the DW data in the LGN ROI using only the original
% tracks.
%
% Franco (C) 2012 Stanford VISTA team.

if notDefined('dtiDataType'); dtiDataType = 'b1000';end
if notDefined('lambda'); lambda = 0;  end
if notDefined('algo');   algo = 'STT';end % TEND or STT
if notDefined('this_roi');   this_roi = 'mct_roi_RTO_12_try2.mat';end% 'mct_r_lgn.mat';end
if notDefined('fitType');   fitType = 'sgd';end
  
saveDir      = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/To12_test';

dataRootPath = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
subfolders   = fullfile('150dirs_b1000_1');
dataDir      = fullfile(dataRootPath,subfolders);
dtFile       = fullfile(dataDir,'dt6.mat');
dwiFile      = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');

% Precomputed ROIs and Fiber groups
roiDir       = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/To12_test';
fiberDir     = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/To12_test';
fgTO12       = fullfile(fiberDir,sprintf('mct_fg_RTO_12_%s.mat',algo));
fgNoILF      = fullfile(fiberDir,sprintf('mct_fg_RTO_12_NO_ILF_%s.mat',algo));

% load the dwi and dti files        
dwi = dwiLoad(dwiFile);
[dtiF, dtiH] = mrDiffusion('off',dtFile);
 
% get te xForms, from a dn to Acpc
xForm.img2acpc = dtiGet(dtiH,'img 2 acpc xform');
xForm.acpc2img = dtiGet(dtiH,'acpc 2 img xform');
close(dtiF), drawnow
clear dtiH dtiF

% load a t1
disp('Loading T1')
nifti = readFileNifti(fullfile(dataRootPath,'t1/t1.nii.gz'));

% Get the ROI
roiName  = fullfile(roiDir,this_roi);
roi      = dtiReadRoi(roiName);
coords = unique(floor(mrAnatXformCoords(xForm.acpc2img ,roi.coords)),'rows');
fprintf('\n[%s] ROI: %s, size(%i,%i)\n',mfilename,roi.name,size(roi.coords))

% Load the contrack fibers and the STT or TEND fibers and add them together.
fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,fgNoILF)
fgNoILF = dtiLoadFiberGroup(fgNoILF);
fgNoILF = dtiXformFiberCoords(fgNoILF, xForm.acpc2img,'img');

fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,fgTO12)
fgTO12 = dtiLoadFiberGroup(fgTO12);
fgTO12  = dtiXformFiberCoords(fgTO12, xForm.acpc2img,'img');

% Clip the fibers to the portion inside the ROI
%fgTO12Clip = mctFGclip(fgTO12,coords);

% Select only the fibers' coordinates going throguht the ROI
% and only the ROI coordinates where a fiber goes through
%unique_fg_coords = fefgGet(fgTO12,'unique image coords');
%coords = intersect(coords,unique_fg_coords,'rows');
fprintf('\n[%s] ROI LGN, size(%i,%i)\n',mfilename,size(coords))


%% (1) Compute the quality fo fit of the combined fiber group.
% Build the MicroTrack model
[Afiber, Aiso, dSig, dSig_demeaned, ~, usedFibers, usedVoxels] = mctBuildDiffusionModel(dwi,fgTO12,coords,[],'ones');
nFiber = size(Afiber,2); 

% build the full model
Afull = [Afiber, Aiso];

% Fit the model.
fiber_wTO12 = mctFitDiffusionModel(Afiber, dSig_demeaned, fitType,lambda);

% Fit the isotropic
iso_w = Aiso \ dSig;
full_w = [fiber_wTO12; iso_w];

% Predict the signal 
pSig_full  = mctComputePredictedSignal(Afull,full_w);

% Compute prediction quality of LiFE WM model
nVoxels = size(Aiso,2); 
nBvecs  = size(Afiber,1) / nVoxels;

% Overall quality across voxels
[~, r2_full] = mctComputePredictionQuality(dSig, pSig_full,2);

% Plot the results
% parse fiber weights and isotropic weights
[wFibers_fu wIsotropic] = mctSortModelWeights(full_w,nFiber);

% Compute the diffusion signal predicted by the original fiber model.
Aorig     = [sum(Afiber,2)  Aiso * iso_w]; 
orig_w    = Aorig \ dSig;
pSig_orig = mctComputePredictedSignal(Aorig,orig_w);
[~, r2_orig] = mctComputePredictionQuality(dSig, pSig_orig,2);

% plot quality of fit.
hQualTO12 = mctPlotQualityOfFit(r2_full,r2_orig, pSig_full,pSig_orig, dSig);

% plotting fiber predictions against the signal  
hFitTO12 = mctDisplayModelFit(wFibers_fu,wIsotropic, dSig, pSig_full, r2_full,  fitType);

% Select the fibers explaining most of the varicance
hits_fibers = fiber_wTO12 >= 0.015; % the fibers with large weights
fa_fibers  = ~hits_fibers;

% now merge again the fibers but do not clip them, We want to show the
% fibers in their full length, to see where they go.
fprintf('\n[%s] Extracting good and bad fibers...\n',mfilename)
fgGood = fgExtract(fgTO12,find(hits_fibers),'keep');
fgGood = dtiXformFiberCoords(fgGood, xForm.img2acpc,'acpc');

fgBad  = fgExtract(fgTO12,find(fa_fibers),'keep');
fgBad  = dtiXformFiberCoords(fgBad, xForm.img2acpc,'acpc');

% The following lines are commented out hoping to not used them for the
% test.
% %% (2) Now repeat the same with fibers not containing the ILF
% % Clip the fibers to the portion inside the ROI
% fgNoILFClip = mctFGclip(fgNoILF,coords);
% 
% %unique_fg_coords = fefgGet(fgNoILF,'unique image coords');
% %coords = intersect(coords,unique_fg_coords,'rows');
% fprintf('\n[%s] ROI LGN size(%i,%i)\n',mfilename,size(coords))
% 
% % Build the MicroTrack model
% [Afiber, Aiso, dSig, dSig_demeaned] = mctBuildDiffusionModel(dwi,fgNoILFClip,coords,[],'ones');
% nFiber = size(Afiber,2); 
% 
% % build the full model
% Afull = [Afiber, Aiso];
% 
% % Fit the model.
% % Find the smallest number fo fibers that explain most of the variance in
% % the diffusion data
% % Fit the fibers
% fiber_wNoILF = mctFitDiffusionModel(Afiber, dSig_demeaned, fitType,lambda);
% 
% % Fit the isotropic
% iso_w = Aiso \ dSig;
% full_w = [fiber_wNoILF; iso_w];
% 
% % Predict the signal 
% pSig_full  = mctComputePredictedSignal(Afull,full_w);
% 
% % Compute prediction quality of LiFE WM model
% nVoxels = size(Aiso,2); 
% nBvecs  = size(Afiber,1) / nVoxels;
% 
% % Compute the diffusion signal predicted by the original fiber model.
% Aorig     = [sum(Afiber,2)  Aiso * iso_w]; 
% orig_w    = Aorig \ dSig;
% pSig_orig = mctComputePredictedSignal(Aorig,orig_w);
% [~, r2_orig] = mctComputePredictionQuality(dSig, pSig_orig,2);
% 
% % Overall quality across voxels
% [~, r2_full] = mctComputePredictionQuality(dSig, pSig_full,2);
% 
% % Plot the results
% % parse fiber weights and isotropic weights
% [wFibers_fu wIsotropic] = mctSortModelWeights(full_w,nFiber);
% 
% % plotting fiber predictions against the signal  
% hFitNoILF = mctDisplayModelFit(wFibers_fu,wIsotropic, dSig, pSig_full, r2_full,  fitType);
% 
% % plot quality of fit.
% hQualNoILF = mctPlotQualityOfFit(r2_full,r2_orig, pSig_full,pSig_orig, dSig);
% 
% % Select the fibers explaining most of the varicance
% hits_fibers = fiber_wNoILF >= 0.015; % the fibers with large weights
% fa_fibers  = ~hits_fibers;
% 
% fprintf('\n[%s] Extracting good and bad fibers...\n',mfilename)
% fgGoodNoILF = fgExtract(fgNoILF,find(hits_fibers),'keep');
% fgBadNoILF  = fgExtract(fgNoILF,find(fa_fibers),'keep');

% Plot fibers with a slice.
h2 = mctNfgDisplayStrands(fgGood,[.75 .2 .3],[],[],[],20,.1.*ones(size(fgGood.fibers)));
hold on
hm = mctDisplayBrainSlice(nifti,[0 0 -7]);
title('Good fibers')

h3 = mctNfgDisplayStrands(fgBad,[.3 .2 .75],[],[],[],20,.1.*ones(size(fgBad.fibers)));
hold on
hm = mctDisplayBrainSlice(nifti,[0 0 -7]);
title('Bad fibers')

h1 = mctNfgDisplayStrands(fgTO12,[.2 .75 .3],[],[],[],20,fiber_wTO12);
hold on
hm = mctDisplayBrainSlice(nifti,[0 0 -7]);
title('All fibers')

keyboard


% Show the fibers
h4 = mctNfgDisplayStrands(fgNoILF,[.2 .75 .3],[],[],[],20,fiber_wNoILF(1:length(fgNoILF.fibers)) .* .1);
h5 = mctNfgDisplayStrands(fgGoodNoILF,[.75 .2 .3],[],[],[],20,.001.*ones(size(fgGoodNoILF.fibers)));
h6 = mctNfgDisplayStrands(fgBadNoILF,[.3 .2 .75],[],[],[],20,.001.*ones(size(fgBadNoILF.fibers)));

fileName  = sprintf('fe_test_TO12_ILF_%s_%s_%s',fitType,dtiDataType,algo);
saveMe = fullfile(saveDir,fileName);
fprintf('\n\nSaving file: %s\n\n',saveMe)
save(saveMe,'-v7.3')

saveas(hFitTO12,fullfile(saveDir,[fileName,'_Fit_TO12_',algo]),'fig')
savefigvista(hFitTO12,[fileName,'_Fit_TO12_',algo],'eps',saveDir,1,1);

saveas(hQualTO12,fullfile(saveDir,[fileName,'_Qual_TO12_',algo]),'fig')
savefigvista(hQualTO12,[fileName,'_Qual_TO12_',algo],'eps',saveDir,1,1);

saveas(hFitNoILF,fullfile(saveDir,[fileName,'_Fit_noILF_',algo]),'fig')
savefigvista(hFitNoILF,[fileName,'_Fit_noILF_',algo],'eps',saveDir,1,1);

saveas(hQualNoILF,fullfile(saveDir,[fileName,'_Qual_noILF_',algo]),'fig')
savefigvista(hQualNoILF,[fileName,'_Qual_noILF_',algo],'eps',saveDir,1,1);

saveas(h1,fullfile(saveDir,[fileName,'_fgAll_TO12_',algo]),'fig')
saveas(h2,fullfile(saveDir,[fileName,'_fgGood_TO12_',algo]),'fig')
saveas(h3,fullfile(saveDir,[fileName,'_fgBad_TO12_',algo]),'fig')

saveas(h4,fullfile(saveDir,[fileName,'_fgAll_noILF_',algo]),'fig')
saveas(h5,fullfile(saveDir,[fileName,'_fgGood_noILF_',algo]),'fig')
saveas(h6,fullfile(saveDir,[fileName,'_fgBad_noILF_',algo]),'fig')


keyboard

end


