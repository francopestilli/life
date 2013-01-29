function s_life_miss_lgn(fitType, lambda,algo, this_roi)
% function s_life_miss_lgn
%
% Load the LGN roi, load some tracks produced with FACT or STT out of such
% roi. 
% 
% Then fits three times to the DW data in the LGN ROI and fit to the
% residual signal left after each fit.
%
% Franco (C) 2012 Stanford VISTA team.

if notDefined('dtiDataType'); dtiDataType = 'b1000';end
if notDefined('lambda'); lambda = 0;end
if notDefined('algo');   algo = 3;end
if notDefined('this_roi');   this_roi = 'mct_r_lgn.mat';end
if notDefined('fitType');   fitType = 'sgd';end

thisDir = pwd;
cd /biac2/wandell6/data/arokem/ModelFits/FP20120420/;
saveDir      = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/LGN_test';
dataRootPath = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
subfolders   = fullfile('150dirs_b1000_1');
dataDir      = fullfile(dataRootPath,subfolders);
dtFile       = fullfile(dataDir,'dt6.mat');
dwiFile      = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');

% PRecomputed ROIs and Fiber groups
roiDir       = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/LGN_test';
fiberDir     = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/LGN_test';
fgLGN        = fullfile(fiberDir,'mct_fg_r_lgn_FACT.mat');

% load the dwi and dti files        
dwi = dwiLoad(dwiFile);
[dtiF, dtiH] = mrDiffusion('off',dtFile);drawnow
 
% get te xForms, from a dn to Acpc
xForm.img2acpc = dtiGet(dtiH,'img 2 acpc xform');
xForm.acpc2img = dtiGet(dtiH,'acpc 2 img xform');
close(dtiF)
clear dtiH dtiF

% Get the ROI
roiName  = fullfile(roiDir,this_roi);
roi      = dtiReadRoi(roiName);
coords = unique(floor(mrAnatXformCoords(xForm.acpc2img ,roi.coords)),'rows');
fprintf('\n[%s] ROI: %s, size(%i,%i)\n',mfilename,roi.name,size(roi.coords))

fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,fgLGN)
fgLgn = dtiLoadFiberGroup(fgLGN);
fgImg = dtiXformFiberCoords(fgLgn, xForm.acpc2img,'img');

% Select only the fibers' coordinates going throguht the ROI
% and only the ROI coordinates where a fiber goes through
unique_fg_coords = fefgGet(fgImg,'unique image coords');
coords = intersect(coords,unique_fg_coords,'rows');
fprintf('\n[%s] ROI IMG %s, size(%i,%i)\n',mfilename,roi.name,size(coords))
clear fg

% Build the MicroTrack model
[Afiber, Aiso, dSig, dSig_demeaned, ~, used_fibers, usedVox] = mctBuildDiffusionModel(dwi,fgImg,coords,[],'ones');
nFiber = size(Afiber,2); 

% build the full model
Afull = [Afiber, Aiso];

% Fit the model.
fiber_w = mctFitDiffusionModel(Afiber, dSig_demeaned, fitType,lambda);
fiber_w = fiber_w(:);

% Fit the isotropic
iso_w = Aiso \ dSig;
full_w = [fiber_w; iso_w];

% Predict the signal 
pSig_full  = mctComputePredictedSignal(Afull,full_w);
pSig_fiber = mctComputePredictedSignal(Afiber,fiber_w);

% Compute prediction quality of LiFE WM model
coords = coords(usedVox,:);
nVoxels = size(Aiso,2); 
nBvecs  = size(Afiber,1) / nVoxels;

% Overall quality across voxels
[rmse_full, r2_full] = mctComputePredictionQuality(dSig, pSig_full);

% Overall quality across voxels
[rmse_fiber, r2_fiber] = mctComputePredictionQuality(dSig_demeaned,pSig_fiber);


% Compute prediction quality of the original WM model
% Compute the diffusion signal predicted by the original fiber model.
Aorig     = [sum(Afiber,2)  Aiso * iso_w]; 
orig_w    = Aorig \ dSig;
pSig_orig = mctComputePredictedSignal(Aorig,orig_w);

clear Afull Afiber Aorig;

% Overall quality across voxels
[rmse_orig, r2_orig] = mctComputePredictionQuality(dSig, pSig_orig);

% parse fiber weights and isotropic weights
[wFibers_fu wIsotropic] = mctSortModelWeights(full_w,nFiber);


%% +++++++++++++++++++++++++++++++++++++++++++++++++++ %%
% Compute residual
fprintf('\n[%s] Working on the second refinement iteration.\n',mfilename);
res = dSig - pSig_full;

% This signal was not explained by the model 
resSig = res + (dSig - dSig_demeaned);

% Reshape the signal to the volume
resSig_img = mctReshape(resSig, nBvecs, nVoxels);

% Load the dwi data.
rawData = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/raw/0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz';
dwi_new = readFileNifti(rawData);

% Add the residuals here into dwRawAligned
% dwiRawAligned can be the original dwiRawAligned lodaed above with the changes in signal at the location of the ROI.
% This is going to be the new dwi that we use for fitting tensors
dwiIter1 = dwi;
for ic = 1:size(coords,1)
  xyz = coords(ic,:);
  % THis is the local dwi file i use.
  dwiIter1.nifti.data(xyz(1),xyz(2),xyz(3),11:end) = resSig_img(:,ic);
  % this si the same data in a format for dtiRawTenorMex.m
  dwi_new.data(xyz(1),xyz(2),xyz(3),11:end)   = resSig_img(:,ic);
end

% Recompute the tensors with the new data. This operation will generate a new dt6 file.
outBaseDir = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/150dirs_b1000_1_iter1/';
dtFileNew = dtiRawFitTensorMex(dwi_new, ...
                                 'raw/0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.bvecs', ...
                                 'raw/0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.bvals', ...
                                 outBaseDir,[],[], [],[],[],1);

% load the new dt6
dt = dtiLoadDt6(dtFileNew);

% Track fibers inside the ROI using an implementation fo FACT.
fgName = sprintf('fe_lgn_fg_iter1_%s',this_roi);
numFibersPerCoord = 1;
fgIter1 = mctTrackInRoi(dt,roi,numFibersPerCoord, fgName,algo);
fgIter1.seeds = [];
fgImgIter1  = dtiXformFiberCoords(fgIter1, xForm.acpc2img,'img');

% Select only the fibers' coordinates going throguht the ROI
% and only the ROI coordinates where a fiber goes through
unique_fg_coordsIter1 = fefgGet(fgImgIter1,'unique image coords');
[coordsIter1, ~, ~] = intersect(coords,unique_fg_coordsIter1,'rows');
fprintf('\n[%s] ROI IMG iter1 %s, size(%i,%i)\n',mfilename,roi.name,size(coordsIter1))

% Build the MicroTrack model
% Afiber is the portion for the fibers
% Aiso is the isotropic portion
% dSig is the diffusion signal (raw)
% dSig_demeaned is the diffusion signal with the mean removed.
[AfiberIter1, AisoIter1, dSigIter1, dSig_demeanedIter1,~,~, usedVoxels] = mctBuildDiffusionModel(dwiIter1,fgImgIter1,coordsIter1,[],'ones');
nFiberIter1 = size(AfiberIter1,2); 

% build the full model
AfullIter1 = [AfiberIter1, AisoIter1];


% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data
% Fit the fibers
fiber_wIter1 = mctFitDiffusionModel(AfiberIter1, dSig_demeanedIter1, fitType,lambda);
fiber_wIter1 = fiber_wIter1(:);

% Fit the isotropic
iso_wIter1 = AisoIter1 \ dSigIter1;
full_wIter1 = [fiber_wIter1; iso_wIter1];

% Predict the signal 
pSig_fullIter1  = mctComputePredictedSignal(AfullIter1,full_wIter1);
pSig_fiberIter1 = mctComputePredictedSignal(AfiberIter1,fiber_wIter1);

% Compute the diffusion signal predicted by the original fiber model.
AorigIter1     = [sum(AfiberIter1,2)  AisoIter1 * iso_wIter1]; 
orig_wIter1    = AorigIter1 \ dSigIter1;
pSig_origIter1 = mctComputePredictedSignal(AorigIter1,orig_wIter1);

% Compute prediction quality of LiFE WM model
% Overall quality across voxels
[~, r2_fullIter1] = mctComputePredictionQuality(dSigIter1, pSig_fullIter1);

% Plot the results
% parse fiber weights and isotropic weights
[wFibers_fuIter1 wIsotropicIter1] = mctSortModelWeights(full_wIter1,nFiberIter1);
wFibers_fiIter1 = mctSortModelWeights(fiber_wIter1,nFiberIter1);

clear AfiberIter1 dwi_new fgIter1 dwiIter1

% plot quality of fit.
% There is good chance the dSig and pSig_orig have different voxels inside.
% So here we clip out the unquated voxels from the signals before computing
% the quality of fit.
% get all the voxels indexes
all_voxIndex  = zeros(size(coords,1),1);
all_voxIndex(usedVoxels) = 1;
all_sigIndex = repmat(all_voxIndex,1,nBvecs)';
all_sigIndex = logical(all_sigIndex(:));


% (2) Test that the sum of the fibers predicts the sum of the R2 from the first fit and the first iteration.
fgImgAll_iter1  = fgMerge(fgImg,fgImgIter1);

[AfiberAll_iter1, AisoAll_iter1, dSigAll_iter1, dSig_demeanedAll_iter1 ~, ~, ~] = mctBuildDiffusionModel(dwi,fgImgAll_iter1,coords,[],'ones');
nFiberAll_iter1 = size(AfiberAll_iter1,2); 

% build the full model
AfullAll_iter1 = [AfiberAll_iter1, AisoAll_iter1];

% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data
% Fit the fibers
fiber_wAll_iter1 = mctFitDiffusionModel(AfiberAll_iter1, dSig_demeanedAll_iter1, fitType,lambda);
fiber_wAll_iter1 = fiber_wAll_iter1(:);

% Fit the isotropic
iso_wAll_iter1 = AisoAll_iter1 \ dSigAll_iter1;
full_wAll_iter1 = [fiber_wAll_iter1; iso_wAll_iter1];

% Predict the signal 
pSig_fullAll_iter1  = mctComputePredictedSignal(AfullAll_iter1,full_wAll_iter1);
pSig_fiberAll_iter1 = mctComputePredictedSignal(AfiberAll_iter1,fiber_wAll_iter1);

%  ##CC## Clear memeory
clear AfullAll_iter1 AfiberAll_iter1 dwi_new fgIter1 dwiIter1

% Compute prediction quality of LiFE WM model
% OverAll_iter1 quality across voxels
[~, r2_fullAll_iter1] = mctComputePredictionQuality(dSigAll_iter1, pSig_fullAll_iter1);

% Plot the results
% parse fiber weights and isotropic weights
[wFibers_fuAll_iter1 wIsotropicAll_iter1] = mctSortModelWeights(full_wAll_iter1,nFiberAll_iter1);



%% ITERATION 2
% Compute residual
fprintf('\n[%s] Working on the second refinement iteration.\n',mfilename);
nVoxels = size(AisoIter1,2); clear AisotIter1

res = [];resSig = []; resSig_img = [];
res = dSigIter1 - pSig_fullIter1;

% This signal was not explained by the model 
resSig = res + (dSigIter1 - dSig_demeanedIter1);

% Reshape the signal to the volume
resSig_img = mctReshape(resSig, nBvecs, nVoxels);

% Load the dwi data.
rawData = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/raw/0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz';
dwi_new = readFileNifti(rawData);

% Add the residuals here into dwRawAligned
% dwiRawAligned can be the original dwiRawAligned lodaed above with the changes in signal at the location of the ROI.
% This is going to be the new dwi that we use for fitting tensors
dwiIter2 = dwi;
for ic = 1:size(coordsIter1,1)
  xyz = coordsIter1(ic,:);
  % THis is the local dwi file i use.
  dwiIter2.nifti.data(xyz(1),xyz(2),xyz(3),11:end) = resSig_img(:,ic);
  % this si the same data in a format for dtiRawTenorMex.m
  dwi_new.data(xyz(1),xyz(2),xyz(3),11:end)   = resSig_img(:,ic);
end

% Recompute the tensors with the new data. This operation will generate a new dt6 file.
outBaseDir = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/150dirs_b1000_1_iter2/';
dtFileNew = dtiRawFitTensorMex(dwi_new, ...
                                 'raw/0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.bvecs', ...
                                 'raw/0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.bvals', ...
                                 outBaseDir,[],[], [],[],[],1);

% load the new dt6
dt = dtiLoadDt6(dtFileNew);

% Track fibers inside the ROI using an implementation fo FACT.
fgName = sprintf('fe_lgn_fg_iter2_%s',this_roi);
numFibersPerCoord = 1;
fgIter2 = mctTrackInRoi(dt,roi,numFibersPerCoord, fgName,algo);
fgIter2.seeds = [];
fgImgIter2  = dtiXformFiberCoords(fgIter2, xForm.acpc2img,'img');


% Select only the fibers' coordinates going throguht the ROI
% and only the ROI coordinates where a fiber goes through
unique_fg_coordsIter2 = fefgGet(fgImgIter2,'unique image coords');
[coordsIter2, ~, ~] = intersect(coordsIter1,unique_fg_coordsIter2,'rows');
% dwiRoi = dwiGet(dwi,'dimage',coords);
fprintf('\n[%s] ROI IMG iter2 %s, size(%i,%i)\n',mfilename,roi.name,size(coordsIter2))

% Build the MicroTrack model
% Afiber is the portion for the fibers
% Aiso is the isotropic portion
% dSig is the diffusion signal (raw)
% dSig_demeaned is the diffusion signal with the mean removed.
[AfiberIter2, AisoIter2, dSigIter2, dSig_demeanedIter2,~,~, usedVoxels] = mctBuildDiffusionModel(dwiIter2,fgImgIter2,coordsIter2,[],'ones');
nFiberIter2 = size(AfiberIter2,2); 

% build the full model
AfullIter2 = [AfiberIter2, AisoIter2];

% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data
% Fit the fibers
fiber_wIter2 = mctFitDiffusionModel(AfiberIter2, dSig_demeanedIter2, fitType,lambda);
fiber_wIter2 = fiber_wIter2(:);

% Fit the isotropic
iso_wIter2 = AisoIter2 \ dSigIter2;
full_wIter2 = [fiber_wIter2; iso_wIter2];

% Predict the signal 
pSig_fullIter2  = mctComputePredictedSignal(AfullIter2,full_wIter2);
pSig_fiberIter2 = mctComputePredictedSignal(AfiberIter2,fiber_wIter2);

% Compute the diffusion signal predicted by the original fiber model.
AorigIter2     = [sum(AfiberIter2,2)  AisoIter2 * iso_wIter2]; 
orig_wIter2    = AorigIter2 \ dSigIter2;
pSig_origIter2 = mctComputePredictedSignal(AorigIter2,orig_wIter2);

% Compute prediction quality of LiFE WM model
% Overall quality across voxels
[~, r2_fullIter2] = mctComputePredictionQuality(dSigIter2, pSig_fullIter2);

% Plot the results
% parse fiber weights and isotropic weights
[wFibers_fuIter2 wIsotropicIter2] = mctSortModelWeights(full_wIter2,nFiberIter2);

clear AfiberIter2 AisoIter2 dwi_new fgIter2 dwiIter2

% plot quality of fit.
% There is good chance the dSig and pSig_orig have different voxels inside.
% So here we clip out the unquated voxels from the signals before computing
% the quality of fit.
% get all the voxels indexes
all_voxIndex  = zeros(size(coordsIter1,1),1);
all_voxIndex(usedVoxels) = 1;
all_sigIndex_iter2 = repmat(all_voxIndex,1,nBvecs)';
all_sigIndex_iter2 = logical(all_sigIndex_iter2(:));

% (2) Test that the sum of the fibers predicts the sum of the R2 from the first fit and the first iteration. 
fgImgAll_iter2  = fgMerge(fgImgAll_iter1,fgImgIter2);

[AfiberAll_iter2, AisoAll_iter2, dSigAll_iter2, dSig_demeanedAll_iter2 ~, ~, ~] = mctBuildDiffusionModel(dwi,fgImgAll_iter2,coords,[],'ones');
nFiberAll_iter2 = size(AfiberAll_iter2,2); 

% build the full model
AfullAll_iter2 = [AfiberAll_iter2, AisoAll_iter2];

% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data
% Fit the fibers
fiber_wAll_iter2 = mctFitDiffusionModel(AfiberAll_iter2, dSig_demeanedAll_iter2, fitType,lambda);
fiber_wAll_iter2 = fiber_wAll_iter2(:);

% Fit the isotropic
iso_wAll_iter2 = AisoAll_iter2 \ dSigAll_iter2;
full_wAll_iter2 = [fiber_wAll_iter2; iso_wAll_iter2];

% Predict the signal 
pSig_fullAll_iter2  = mctComputePredictedSignal(AfullAll_iter2,full_wAll_iter2);

%  ##CC## Clear memeory
clear AfullAll_iter2 AfiberAll_iter2 dwi_new fgIter1 dwiIter1

% Compute prediction quality of LiFE WM model
% OverAll_iter2 quality across voxels
[~, r2_fullAll_iter2] = mctComputePredictionQuality(dSigAll_iter2, pSig_fullAll_iter2);

% Plot the results
% parse fiber weights and isotropic weights
[wFibers_fuAll_iter2 wIsotropicAll_iter2] = mctSortModelWeights(full_wAll_iter2,nFiberAll_iter2);



%% Select the fibers explaining most of the varicance
hits_fibers = fiber_w >= 0.001; % the fibers with large weights
fa_fibers  = ~hits_fibers;

fgGood = fgExtract(fgImg,find(hits_fibers),'keep');
fgGood.seeds = [];
fgBad  = fgExtract(fgImg,find(fa_fibers),'keep');
faBad.seeds = [];

% Now clip the original and the new fibers tp shwo the local difference in
% structure evidenced by LiFE
fgGoodClip = mctFGclip(fgGood,coords);
fgImg.seeds = [];
fgImgClip = mctFGclip(fgImg,coords);


% save results to disk
fileName  = sprintf('%s_%s_%s_lambda%i_algo%i',this_roi(1:end-4),fitType,dtiDataType,10000*lambda,algo);
saveMe = fullfile(saveDir,fileName);
fprintf('\n\nSaving file: %s\n\n',saveMe)
save(saveMe,'-v7.3')

% Save figures
fileName  = sprintf('%s_%s_%s_lambda%i_algo%i',this_roi(1:end-4),fitType,dtiDataType,lambda,algo);

h18 = mctNfgDisplayStrands(fgGoodClip,[.5 .75 .3],[],[],[],20,.00001.*ones(size(fgGoodClip.fibers)));
h19 = mctNfgDisplayStrands(fgImgClip,[.75 .5 .3],[],[],[],20,.00001.*ones(size(fgImgClip.fibers)));

% plotting fiber predictions against the signal  
h1  = mctDisplayModelFit(wFibers_fu,wIsotropic, dSig, pSig_full, r2_full,  fitType);
h2  = mctPlotQualityOfFit(r2_full,r2_orig, pSig_full,pSig_orig, dSig);
h3  = mctNfgDisplayStrands(fgImg,[.5 .75 .3],[],[],[],20,.001.*ones(size(fgImg.fibers)));

h4  = mctDisplayModelFit(wFibers_fuIter1,wIsotropicIter1, dSigIter1, pSig_fullIter1, r2_fullIter1,  fitType);
h5  = mctPlotQualityOfFit(r2_fullIter1,r2_orig, pSig_fullIter1,pSig_origIter1, resSig(all_sigIndex));
h6  = mctNfgDisplayStrands(fgImgIter1,[.5 .75 .3],[],[],[],20,.001.*ones(size(fgImgIter1.fibers)));

h7  = mctDisplayModelFit(wFibers_fuIter2,wIsotropicIter2, dSigIter2, pSig_fullIter2, r2_fullIter2,  fitType);
h8  = mctPlotQualityOfFit(r2_fullIter2,r2_orig, pSig_fullIter2,pSig_origIter2, resSig(all_sigIndex_iter2));
h9  = mctNfgDisplayStrands(fgImgIter2,[.5 .75 .3],[],[],[],20,.001.*ones(size(fgImgIter2.fibers)));

h10  = mctDisplayModelFit(wFibers_fuAll_iter1,wIsotropicAll_iter1, dSigAll_iter1, pSig_fullAll_iter1, r2_fullAll_iter1,  fitType);
h11  = mctPlotQualityOfFit(r2_fullAll_iter1,r2_orig, pSig_fullAll_iter1,pSig_orig, dSigAll_iter1);
h12  = mctNfgDisplayStrands(fgImgAll_iter1,[.5 .75 .3],[],[],[],20,.001.*ones(size(fgImgAll_iter1.fibers)));

h13  = mctDisplayModelFit(wFibers_fuAll_iter2,wIsotropicAll_iter2, dSigAll_iter2, pSig_fullAll_iter2, r2_fullAll_iter2,  fitType);
h14  = mctPlotQualityOfFit(r2_fullAll_iter2,r2_orig, pSig_fullAll_iter2,pSig_orig, dSigAll_iter2);
h15  = mctNfgDisplayStrands(fgImgAll_iter2,[.5 .75 .3],[],[],[],20,.001.*ones(size(fgImgAll_iter2.fibers)));

h16 = mctNfgDisplayStrands(fgGood,[.5 .75 .3],[],[],[],20,.001.*ones(size(fgGood.fibers)));
h17 = mctNfgDisplayStrands(fgBad,[.75 .5 .3],[],[],[],20,.001.*ones(size(fgBad.fibers)));


% Save figures
savefigvista(h1,[fileName,'_weigths_1'],'eps',saveDir,'/',1,1);
savefigvista(h2,[fileName,'_r2_1'],'eps',saveDir,'/',1,1);
saveas(      h3,fullfile(saveDir,[fileName,'_fg_1']),'fig')

savefigvista(h4,[fileName,'_weigths_2'],'eps',saveDir,'/',1,1);
savefigvista(h5,[fileName,'_r2_2'],'eps',saveDir,'/',1,1);
saveas(      h6,fullfile(saveDir,[fileName,'_fg_2']),'fig')

savefigvista(h7,[fileName,'_weigths_3'],'eps',saveDir,'/',1,1);
savefigvista(h8,[fileName,'_r2_3'],'eps',saveDir,'/',1,1);
saveas(      h9,fullfile(saveDir,[fileName,'_fg_3']),'fig')

savefigvista(h10,[fileName,'_weigths_all1'],'eps',saveDir,'/',1,1);
savefigvista(h11,[fileName,'_r2_all1'],'eps',saveDir,'/',1,1);
saveas(      h12,fullfile(saveDir,[fileName,'_fg_all1']),'fig')

savefigvista(h13,[fileName,'_weigths_all2'],'eps',saveDir,'/',1,1);
savefigvista(h14,[fileName,'_r2_all2'],'eps',saveDir,'/',1,1);
saveas(      h15,fullfile(saveDir,[fileName,'_fg_all2']),'fig')

saveas(      h16,fullfile(saveDir,[fileName,'_fg_good']),'fig')
saveas(      h17,fullfile(saveDir,[fileName,'_fg_bad']),'fig')

saveas(      h18,fullfile(saveDir,[fileName,'_fg_good_clipped']),'fig')
saveas(      h19,fullfile(saveDir,[fileName,'_fg_orig_clipped']),'fig')

keyboard
cd(thisDir)

