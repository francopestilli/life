function s_life_miss(fitType, lambda,algo, this_roi)
% function s_life_miss
%
% This script loads two rois, track with fact and fits both of them with microtrack.
%
% - The corpus callosum ROI should have a good fit with microtrack, because
% all the fibers go in the same direction.
%
% - The prefrontal ROi should have a worse fit, because there will be
% fibers going in different directions, directions not accounted for bu the
% deterministic algorithm of FACT.
% 
% Franco
%       
% FP rois: 
% 'mct_callosal_1mm.mat'; 'mct_callosal_4mm.mat';
% 'mct_callosal_2mm.mat'; 'mct_prefrontal_1mm.mat';
% 'mct_prefrontal_10mm.mat'; 'mct_prefrontal_4mm.mat';
% 'mct_callosal_40mm.mat'; 'mct_visual_20mm.mat'; 'mct_prefrontal_2mm.mat';
% 'mct_rt12.mat';'mct_or.mat';'mct_occipital_lobe.mat';
%
% (C) 2012 Stanford VISTA team.

if notDefined('dtiDataType'); dtiDataType = 'b1000';end
if notDefined('lambda'); lambda = 0;end
if notDefined('algo');   algo = 3;end
if notDefined('this_roi');   this_roi = 'mct_prefrontal_10mm.mat';end
if notDefined('fitType');   fitType = 'sgd';end

dataRootPath = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
subfolders   = fullfile('150dirs_b1000_1');
dataDir      = fullfile(dataRootPath,subfolders);
roiDir       = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
dtFile       = fullfile(dataDir,'dt6.mat');
dwiFile      = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');

cd /biac2/wandell6/data/arokem/ModelFits/FP20120420/

% load the dwi and dti files        
dwi = dwiLoad(dwiFile);

% dti6
dt           = dtiLoadDt6(dtFile);

[dtiF, dtiH] = mrDiffusion('off',dtFile);

% get te xForms, from a dn to Acpc
xForm.img2acpc = dtiGet(dtiH,'img 2 acpc xform');
xForm.acpc2img = dtiGet(dtiH,'acpc 2 img xform');

%% Get the ROI
roiName  = fullfile(roiDir,'ROIs',this_roi);
roi      = dtiReadRoi(roiName);
fprintf('\n[%s] ROI: %s, size(%i,%i)\n',mfilename,roi.name,size(roi.coords))

% add the roi to the dtiH
dtiH = dtiSet(dtiH,'add roi',roi);

%% Track fibers inside the ROI using an implementation fo FACT.
fgName = sprintf('%s', '/fibers/WholeBrainFG.mat');

if ~(exist(fullfile(dataDir,fgName),'file') == 2)
    numFibersPerCoord = 1;
    % 0=FACT Euler, 1=FACT RK4, 2=TEND Euler, 3=TEND RK4,
    fg = mctTrackInRoi(dt,roi,numFibersPerCoord, fgName,algo);
    fg.seeds = [];
    
    % write the fiber group to file    
    %fprintf('\n[%s]\n saving fibers to file: %s\n',mfilename,fullfile(dataDir,fgName))
    %dtiWriteFiberGroup(fg,fullfile(dataDir,fgName))
    %mtrExportFibers(fg,'fgName',xForm.img2acpc);

elseif exist(fullfile(dataDir,fgName),'file')
    fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,fullfile(dataDir,fgName))
    fg = dtiLoadFiberGroup(fullfile(dataDir,fgName));

else
keyboard

end

% ##CC## clear memory.
clear dtiH dtiF dt

%% Transform the fibers and ROI coordinates in image coordinates
coords = unique(floor(mrAnatXformCoords(xForm.acpc2img ,roi.coords)),'rows');
fgImg  = dtiXformFiberCoords(fg, xForm.acpc2img,'img');

% Select only the fibers' coordinates going throguht the ROI
% and only the ROI coordinates where a fiber goes through
unique_fg_coords = fefgGet(fgImg,'unique image coords');
coords = intersect(coords,unique_fg_coords,'rows');
% dwiRoi = dwiGet(dwi,'dimage',coords);
fprintf('\n[%s] ROI IMG %s, size(%i,%i)\n',mfilename,roi.name,size(coords))

% Clip the fibers to the portion inside the ROI
%fgImg = mctFGclip(fgImg,coords);

% ##CC## Clear memmeory
clear fg

%% Build the MicroTrack model
% Afiber is the portion for the fibers
% Aiso is the isotropic portion
% dSig is the diffusion signal (raw)
% dSig_demeaned is the diffusion signal with the mean removed.
[Afiber, Aiso, dSig, dSig_demeaned, ~, used_fibers, usedVox] = mctBuildDiffusionModel(dwi,fgImg,coords,[],'ones');
nFiber = size(Afiber,2); 

keyboard

% build the full model
Afull = [Afiber, Aiso];

%% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data
% Fit the fibers
fiber_w = mctFitDiffusionModel(Afiber, dSig_demeaned, fitType,lambda);
fiber_w = fiber_w(:);

% Fit the isotropic
iso_w = Aiso \ dSig;
full_w = [fiber_w; iso_w];

%% Predict the signal 
pSig_full  = mctComputePredictedSignal(Afull,full_w);
pSig_fiber = mctComputePredictedSignal(Afiber,fiber_w);


%% Compute prediction quality of LiFE WM model
coords = coords(usedVox,:);
nVoxels = size(Aiso,2); 
nBvecs  = size(Afiber,1) / nVoxels;

% Overall quality across voxels
[rmse_full, r2_full] = mctComputePredictionQuality(dSig, pSig_full);

% Overall quality across voxels
[rmse_fiber, r2_fiber] = mctComputePredictionQuality(dSig_demeaned,pSig_fiber);


%% Compute prediction quality of the original WM model
% Compute the diffusion signal predicted by the original fiber model.
Aorig     = [sum(Afiber,2)  Aiso * iso_w]; 
orig_w    = Aorig \ dSig;
pSig_orig = mctComputePredictedSignal(Aorig,orig_w);

% ##CC## clear memory from all the variables i will not use again.
clear Afull Afiber Aorig;

% Overall quality across voxels
[rmse_orig, r2_orig] = mctComputePredictionQuality(dSig, pSig_orig);

% parse fiber weights and isotropic weights
[wFibers_fu wIsotropic] = mctSortModelWeights(full_w,nFiber);
wFibers_fi = mctSortModelWeights(fiber_w,nFiber);

% plotting fiber predictions against the signal  
h1 = mctDisplayModelFit(wFibers_fu,wIsotropic, dSig, pSig_full, r2_full,  fitType);

% plot quality of fit.
h2 = mctPlotQualityOfFit(r2_full,r2_orig, pSig_full,pSig_orig, dSig);

% Plot fibers.
h3 = mctNfgDisplayStrands(fgImg);

% Save figures
saveDir = '/biac2/wandell6/data/frk/LiFE';
fileName  = sprintf('%s_%s_%s_lambda%i_algo%i',this_roi(1:end-4),fitType,dtiDataType,lambda,algo);

savefigvista(h1,[fileName,'_weigths1'],'eps',saveDir,'/',1,1);
savefigvista(h2,[fileName,'_r21'],'eps',saveDir,'/',1,1);
saveas(h3,fullfile(saveDir,[fileName,'_fg1']),'fig')

% +++++++++++++++++++++++++++++++++++++++++++++++++++ %%
%% Compute residual
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

%% Track fibers inside the ROI using an implementation fo FACT.
fgName = sprintf('fact_fibers_iter1_%s',this_roi);
numFibersPerCoord = 1;
fgIter1 = mctTrackInRoi(dt,roi,numFibersPerCoord, fgName,algo);
fgIter1.seeds = [];
fgImgIter1  = dtiXformFiberCoords(fgIter1, xForm.acpc2img,'img');

% Select only the fibers' coordinates going throguht the ROI
% and only the ROI coordinates where a fiber goes through
unique_fg_coordsIter1 = fefgGet(fgImgIter1,'unique image coords');
[coordsIter1, usedVoxels, usedFibers] = intersect(coords,unique_fg_coordsIter1,'rows');
% dwiRoi = dwiGet(dwi,'dimage',coords);
fprintf('\n[%s] ROI IMG iter1 %s, size(%i,%i)\n',mfilename,roi.name,size(coordsIter1))

% Clip the fibers to the portion inside the ROI
fgImgIter1 = mctFGclip(fgImgIter1,coordsIter1);


%% Build the MicroTrack model
% Afiber is the portion for the fibers
% Aiso is the isotropic portion
% dSig is the diffusion signal (raw)
% dSig_demeaned is the diffusion signal with the mean removed.
[AfiberIter1, AisoIter1, dSigIter1, dSig_demeanedIter1,~,usedFibers, usedVoxels] = mctBuildDiffusionModel(dwiIter1,fgImgIter1,coordsIter1,[],'ones');
nFiberIter1 = size(AfiberIter1,2); 

% build the full model
AfullIter1 = [AfiberIter1, AisoIter1];


%% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data
% Fit the fibers
fiber_wIter1 = mctFitDiffusionModel(AfiberIter1, dSig_demeanedIter1, fitType,lambda);
fiber_wIter1 = fiber_wIter1(:);

% Fit the isotropic
iso_wIter1 = AisoIter1 \ dSigIter1;
full_wIter1 = [fiber_wIter1; iso_wIter1];

%% Predict the signal 
pSig_fullIter1  = mctComputePredictedSignal(AfullIter1,full_wIter1);
pSig_fiberIter1 = mctComputePredictedSignal(AfiberIter1,fiber_wIter1);

%% Compute the diffusion signal predicted by the original fiber model.
AorigIter1     = [sum(AfiberIter1,2)  AisoIter1 * iso_wIter1]; 
orig_wIter1    = AorigIter1 \ dSigIter1;
pSig_origIter1 = mctComputePredictedSignal(AorigIter1,orig_wIter1);

%% Compute prediction quality of LiFE WM model
% Overall quality across voxels
[rmse_fullIter1, r2_fullIter1] = mctComputePredictionQuality(dSigIter1, pSig_fullIter1);

% Overall quality across voxels
[rmse_fiberIter1, r2_fiberIter1] = mctComputePredictionQuality(dSig_demeanedIter1,pSig_fiberIter1);

%% Plot the results
% parse fiber weights and isotropic weights
[wFibers_fuIter1 wIsotropicIter1] = mctSortModelWeights(full_wIter1,nFiberIter1);
wFibers_fiIter1 = mctSortModelWeights(fiber_wIter1,nFiberIter1);

% plotting fiber predictions against the signal  
h4 = mctDisplayModelFit(wFibers_fuIter1,wIsotropicIter1, dSigIter1, pSig_fullIter1, r2_fullIter1,  fitType);

%  ##CC## Clear memeory
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

h5 = mctPlotQualityOfFit(r2_fullIter1,r2_orig, pSig_fullIter1,pSig_origIter1, resSig(all_sigIndex));

% Plot fibers.
h6 = mctNfgDisplayStrands(fgImgIter1);

% Save figures
saveDir = '/biac2/wandell6/data/frk/LiFE';
fileName  = sprintf('%s_%s_%s_lambda%i_algo%i',this_roi(1:end-4),fitType,dtiDataType,lambda,algo);

savefigvista(h4,[fileName,'_weigths2'],'eps',saveDir,'/',1,1);
savefigvista(h5,[fileName,'_r22'],'eps',saveDir,'/',1,1);
saveas(h6,fullfile(saveDir,[fileName,'_fg2']),'fig')

%_________________________________________________________________________%
%% (2) Test that the sum of the fibers predicts the sum of the R2 from the first fit and the first iteration.
fgImgAll  = fgMerge(fgImg,fgImgIter1);
%fgImgAll  = fgMerge(fgImgAll,fgImgIter2);

[AfiberAll, AisoAll, dSigAll, dSig_demeanedAll ~, used_fibersAll, usedVoxAll] = mctBuildDiffusionModel(dwi,fgImgAll,coords,[],'ones');
nFiberAll = size(AfiberAll,2); 

% build the full model
AfullAll = [AfiberAll, AisoAll];

%% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data
% Fit the fibers
fiber_wAll = mctFitDiffusionModel(AfiberAll, dSig_demeanedAll, fitType,lambda);
fiber_wAll = fiber_wAll(:);

% Fit the isotropic
iso_wAll = AisoAll \ dSigAll;
full_wAll = [fiber_wAll; iso_wAll];

%% Predict the signal 
pSig_fullAll  = mctComputePredictedSignal(AfullAll,full_wAll);
pSig_fiberAll = mctComputePredictedSignal(AfiberAll,fiber_wAll);

%  ##CC## Clear memeory
clear AfullAll AfiberAll dwi_new fgIter1 dwiIter1

%% Compute prediction quality of LiFE WM model
% Overall quality across voxels
[rmse_fullAll, r2_fullAll] = mctComputePredictionQuality(dSigAll, pSig_fullAll);

% Overall quality across voxels
[rmse_fiberAll, r2_fiberAll] = mctComputePredictionQuality(dSig_demeanedAll,pSig_fiberAll);

%% Plot the results
% parse fiber weights and isotropic weights
[wFibers_fuAll wIsotropicAll] = mctSortModelWeights(full_wAll,nFiberAll);
wFibers_fiAll = mctSortModelWeights(fiber_wAll,nFiberAll);

% plotting fiber predictions against the signal  
h7 = mctDisplayModelFit(wFibers_fuAll,wIsotropicAll, dSigAll, pSig_fullAll, r2_fullAll,  fitType);

% plot quality of fit.
h8 = mctPlotQualityOfFit(r2_fullAll,r2_orig, pSig_fullAll,pSig_orig, dSigAll);

% Plot fibers.
h9 = mctNfgDisplayStrands(fgImgAll);

% Save figures
savefigvista(h7,[fileName,'_weigths3'],'eps',saveDir,'/',1,1);
savefigvista(h8,[fileName,'_r23'],'eps',saveDir,'/',1,1);
saveas(h9,fullfile(saveDir,[fileName,'_fg3']),'fig')


%% Select the fibers explaining most of the varicance
hits_fibers = fiber_w >= 0.001; % the fibers with large weights
fa_fibers  = ~hits_fibers;

fgGood = fgExtract(fgImg,find(hits_fibers),'keep');
fgBad  = fgExtract(fgImg,find(fa_fibers),'keep');

h10 = mctNfgDisplayStrands(fgGood);
h11 = mctNfgDisplayStrands(fgBad);

% Save figures
saveas(h10,fullfile(saveDir,[fileName,'_fg3good']),'fig')
saveas(h11,fullfile(saveDir,[fileName,'_fg3bad']),'fig')

% save results to disk
saveDir = '/biac2/wandell6/data/frk/LiFE';
fileName  = sprintf('%s_%s_%s_lambda%i_algo%i',this_roi(1:end-4),fitType,dtiDataType,10000*lambda,algo);
saveMe = fullfile(saveDir,fileName);
fprintf('\n\nSaving file: %s\n\n',saveMe)
save(saveMe,'-v7.3')

end


% 
% %%------------------------------------------------------%%
% %% ITERATION 2
% 
% % +++++++++++++++++++++++++++++++++++++++++++++++++++ %%
% %% Compute residual
% fprintf('\n[%s] Working on the second refinement iteration.\n',mfilename);
% nVoxels = size(AisoIter1,2); clear AisotIter1
% 
% res = [];resSig = []; resSig_img = [];
% res = dSigIter1 - pSig_fullIter1;
% 
% % This signal was not explained by the model 
% resSig = res + (dSigIter1 - dSig_demeanedIter1);
% 
% % Reshape the signal to the volume
% resSig_img = mctReshape(resSig, nBvecs, nVoxels);
% 
% % Load the dwi data.
% rawData = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/raw/0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz';
% dwi_new = readFileNifti(rawData);
% 
% % Add the residuals here into dwRawAligned
% % dwiRawAligned can be the original dwiRawAligned lodaed above with the changes in signal at the location of the ROI.
% % This is going to be the new dwi that we use for fitting tensors
% dwiIter2 = dwi;
% for ic = 1:size(coordsIter1,1)
%   xyz = coordsIter1(ic,:);
%   % THis is the local dwi file i use.
%   dwiIter2.nifti.data(xyz(1),xyz(2),xyz(3),11:end) = resSig_img(:,ic);
%   % this si the same data in a format for dtiRawTenorMex.m
%   dwi_new.data(xyz(1),xyz(2),xyz(3),11:end)   = resSig_img(:,ic);
% end
% 
% % Recompute the tensors with the new data. This operation will generate a new dt6 file.
% outBaseDir = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/150dirs_b1000_1_iter2/';
% dtFileNew = dtiRawFitTensorMex(dwi_new, ...
%                                  'raw/0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.bvecs', ...
%                                  'raw/0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.bvals', ...
%                                  outBaseDir,[],[], [],[],[],1);
% 
% 
% 
% % load the new dt6
% dt = dtiLoadDt6(dtFileNew);
% 
% %% Track fibers inside the ROI using an implementation fo FACT.
% fgName = sprintf('fact_fibers_iter2_%s',this_roi);
% numFibersPerCoord = 1;
% fgIter2 = mctTrackInRoi(dt,roi,numFibersPerCoord, fgName,algo);
% fgIter2.seeds = [];
% fgImgIter2  = dtiXformFiberCoords(fgIter2, xForm.acpc2img,'img');
% 
% % Select only the fibers' coordinates going throguht the ROI
% % and only the ROI coordinates where a fiber goes through
% unique_fg_coordsIter2 = fefgGet(fgImgIter2,'unique image coords');
% [coordsIter2, usedVoxels, usedFibers] = intersect(coordsIter1,unique_fg_coordsIter2,'rows');
% % dwiRoi = dwiGet(dwi,'dimage',coords);
% fprintf('\n[%s] ROI IMG iter2 %s, size(%i,%i)\n',mfilename,roi.name,size(coordsIter2))
% 
% % Clip the fibers to the portion inside the ROI
% fgImgIter2 = mctFGclip(fgImgIter2,coordsIter2);
% 
% 
% %% Build the MicroTrack model
% % Afiber is the portion for the fibers
% % Aiso is the isotropic portion
% % dSig is the diffusion signal (raw)
% % dSig_demeaned is the diffusion signal with the mean removed.
% [AfiberIter2, AisoIter2, dSigIter2, dSig_demeanedIter2,~,usedFibers, usedVoxels] = mctBuildDiffusionModel(dwiIter2,fgImgIter2,coordsIter2,[],'ones');
% nFiberIter2 = size(AfiberIter2,2); 
% 
% % build the full model
% AfullIter2 = [AfiberIter2, AisoIter2];
% 
% 
% %% Fit the model.
% % Find the smallest number fo fibers that explain most of the variance in
% % the diffusion data
% % Fit the fibers
% fiber_wIter2 = mctFitDiffusionModel(AfiberIter2, dSig_demeanedIter2, fitType,lambda);
% fiber_wIter2 = fiber_wIter2(:);
% 
% % Fit the isotropic
% iso_wIter2 = AisoIter2 \ dSigIter2;
% full_wIter2 = [fiber_wIter2; iso_wIter2];
% 
% %% Predict the signal 
% pSig_fullIter2  = mctComputePredictedSignal(AfullIter2,full_wIter2);
% pSig_fiberIter2 = mctComputePredictedSignal(AfiberIter2,fiber_wIter2);
% 
% %% Compute the diffusion signal predicted by the original fiber model.
% AorigIter2     = [sum(AfiberIter2,2)  AisoIter2 * iso_wIter2]; 
% orig_wIter2    = AorigIter2 \ dSigIter2;
% pSig_origIter2 = mctComputePredictedSignal(AorigIter2,orig_wIter2);
% 
% %% Compute prediction quality of LiFE WM model
% % Overall quality across voxels
% [rmse_fullIter2, r2_fullIter2] = mctComputePredictionQuality(dSigIter2, pSig_fullIter2);
% 
% % Overall quality across voxels
% [rmse_fiberIter2, r2_fiberIter2] = mctComputePredictionQuality(dSig_demeanedIter2,pSig_fiberIter2);
% 
% %% Plot the results
% % parse fiber weights and isotropic weights
% [wFibers_fuIter2 wIsotropicIter2] = mctSortModelWeights(full_wIter2,nFiberIter2);
% wFibers_fiIter2 = mctSortModelWeights(fiber_wIter2,nFiberIter2);
% 
% % plotting fiber predictions against the signal  
% mctDisplayModelFit(wFibers_fuIter2,wIsotropicIter2, dSigIter2, pSig_fullIter2, r2_fullIter2,  fitType);
% 
% %  ##CC## Clear memeory
% clear AfiberIter2 AisoIter2 dwi_new fgIter2 dwiIter2
% 
% % plot quality of fit.
% % There is good chance the dSig and pSig_orig have different voxels inside.
% % So here we clip out the unquated voxels from the signals before computing
% % the quality of fit.
% % get all the voxels indexes
% all_voxIndex  = zeros(size(coordsIter1,1),1);
% all_voxIndex(usedVoxels) = 1;
% all_sigIndex = repmat(all_voxIndex,1,nBvecs)';
% all_sigIndex = logical(all_sigIndex(:));
% 
% mctPlotQualityOfFit(r2_fullIter2,r2_orig, pSig_fullIter2,pSig_origIter2, resSig(all_sigIndex));
% 
% % Plot fibers.
% mctNfgDisplayStrands(fgImgIter2);
% 
% 
% 