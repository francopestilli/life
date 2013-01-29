function s_life_lgn_vi(fitType, lambda,algo, this_roi)
% function s_life_lgn_vi
%
% Load the LGn roi, load some tracks produced with STT or STT out of such
% roi. 
% 
% Then load tracks of the optic radiation tracked using Contrack.
%
% Then produce a fit to the DW data in the LGN ROI using the original
% tracks summed with the Optic Radiation tracks.
%
% Then produce a fit to the DW data in the LGN ROI using only the original
% tracks.
%
% Franco (C) 2012 Stanford VISTA team.

if notDefined('dtiDataType'); dtiDataType = 'b1000';end
if notDefined('lambda'); lambda = 0;end
if notDefined('algo');   algo = 3;end
if notDefined('this_roi');   this_roi = 'mct_ROI_white_matter_around_lgn.mat';end% 'mct_r_lgn.mat';end
if notDefined('fitType');   fitType = 'sgd';end

saveDir      = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/LGN_test';

dataRootPath = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
subfolders   = fullfile('150dirs_b1000_1');
dataDir      = fullfile(dataRootPath,subfolders);
dtFile       = fullfile(dataDir,'dt6.mat');
dwiFile      = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');

% PRecomputed ROIs and Fiber groups
roiDir       = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/LGN_test';
fiberDir     = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/LGN_test';
fgLGN        = fullfile(fiberDir,'mct_FG_white_matter_around_lgn_STT.mat');
fgContrack = fullfile(fiberDir,'mct_fg_r_optic_radiation_contrack.mat');

% load the dwi and dti files        
dwi = dwiLoad(dwiFile);
dt           = dtiLoadDt6(dtFile);
[dtiF, dtiH] = mrDiffusion('off',dtFile);
 
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

% Load the contrack fibers and the STT or TEND fibers and add them together.
fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,fgContrack)
fgCon = dtiLoadFiberGroup(fgContrack);
fgCon  = dtiXformFiberCoords(fgCon, xForm.acpc2img,'img');

fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,fgLGN)
fgLgn = dtiLoadFiberGroup(fgLGN);
fgLgn  = dtiXformFiberCoords(fgLgn, xForm.acpc2img,'img');

% Merge the two ROIs
fgImg = fgMerge(fgLgn,fgCon,'lgn_opritc_radiation');

% Clip the fibers to the portion inside the ROI
fgImg = mctFGclip(fgImg,coords);

% Select only the fibers' coordinates going throguht the ROI
% and only the ROI coordinates where a fiber goes through
unique_fg_coords = fefgGet(fgImg,'unique image coords');
coords = intersect(coords,unique_fg_coords,'rows');
fprintf('\n[%s] ROI LGN, size(%i,%i)\n',mfilename,size(coords))

%% (1) Compute the quality fo fit of the combined fiber group.
% Build the MicroTrack model
[Afiber, Aiso, dSig, dSig_demeaned] = mctBuildDiffusionModel(dwi,fgImg,coords,[],'ones');
nFiber = size(Afiber,2); 

% build the full model
Afull = [Afiber, Aiso];

% Fit the model.
fiber_w = mctFitDiffusionModel(Afiber, dSig_demeaned, fitType,lambda);

% Fit the isotropic
iso_w = Aiso \ dSig;
full_w = [fiber_w; iso_w];

% Predict the signal 
pSig_full  = mctComputePredictedSignal(Afull,full_w);
pSig_fiber = mctComputePredictedSignal(Afiber,fiber_w);

% Compute prediction quality of LiFE WM model
nVoxels = size(Aiso,2); 
nBvecs  = size(Afiber,1) / nVoxels;

% Overall quality across voxels
[rmse_full, r2_full] = mctComputePredictionQuality(dSig, pSig_full);

% Overall quality across voxels
[rmse_fiber, r2_fiber] = mctComputePredictionQuality(dSig_demeaned,pSig_fiber);

% Plot the results
% parse fiber weights and isotropic weights
[wFibers_fu wIsotropic] = mctSortModelWeights(full_w,nFiber);

% plotting fiber predictions against the signal  
hFitCon = mctDisplayModelFit(wFibers_fu,wIsotropic, dSig, pSig_full, r2_full,  fitType);

% Compute the diffusion signal predicted by the original fiber model.
Aorig     = [sum(Afiber,2)  Aiso * iso_w]; 
orig_w    = Aorig \ dSig;
pSig_orig = mctComputePredictedSignal(Aorig,orig_w);
[rmse_orig, r2_orig] = mctComputePredictionQuality(dSig, pSig_orig);

% plot quality of fit.
hQualCon = mctPlotQualityOfFit(r2_full,r2_orig, pSig_full,pSig_orig, dSig);

% Select the fibers explaining most of the varicance
hits_fibers = fiber_w >= 0.00001; % the fibers with large weights
fa_fibers  = ~hits_fibers;

% now merge again the fibers but do not clip them, We want to show the
% fibers in their full length, to see where they go.
fprintf('\n[%s] Extracting good and bad fibers...\n',mfilename)
fgGoodCon = fgExtract(fgImg,find(hits_fibers(1:length(fgImg.fibers))),'keep');
fgBadCon  = fgExtract(fgImg,find(fa_fibers(1:length(fgImg.fibers))),'keep');


%% (2) Now repeat the same but removing the optic radiation.
% Clip the fibers to the portion inside the ROI
fgLgnImg = mctFGclip(fgLgn,coords);

unique_fg_coords = fefgGet(fgLgnImg,'unique image coords');
coords = intersect(coords,unique_fg_coords,'rows');
fprintf('\n[%s] ROI LGN size(%i,%i)\n',mfilename,size(coords))

% Build the MicroTrack model
[Afiber, Aiso, dSig, dSig_demeaned] = mctBuildDiffusionModel(dwi,fgLgnImg,coords,[],'ones');
nFiber = size(Afiber,2); 

% build the full model
Afull = [Afiber, Aiso];

% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data
% Fit the fibers
fiber_w = mctFitDiffusionModel(Afiber, dSig_demeaned, fitType,lambda);

% Fit the isotropic
iso_w = Aiso \ dSig;
full_w = [fiber_w; iso_w];

% Predict the signal 
pSig_full  = mctComputePredictedSignal(Afull,full_w);
pSig_fiber = mctComputePredictedSignal(Afiber,fiber_w);

% Compute prediction quality of LiFE WM model
nVoxels = size(Aiso,2); 
nBvecs  = size(Afiber,1) / nVoxels;

% Compute the diffusion signal predicted by the original fiber model.
Aorig     = [sum(Afiber,2)  Aiso * iso_w]; 
orig_w    = Aorig \ dSig;
pSig_orig = mctComputePredictedSignal(Aorig,orig_w);
[rmse_orig, r2_orig] = mctComputePredictionQuality(dSig, pSig_orig);

% Overall quality across voxels
[rmse_full, r2_full] = mctComputePredictionQuality(dSig, pSig_full);

% Overall quality across voxels
[rmse_fiber, r2_fiber] = mctComputePredictionQuality(dSig_demeaned,pSig_fiber);

% Plot the results
% parse fiber weights and isotropic weights
[wFibers_fu wIsotropic] = mctSortModelWeights(full_w,nFiber);

% plotting fiber predictions against the signal  
hFitNocon = mctDisplayModelFit(wFibers_fu,wIsotropic, dSig, pSig_full, r2_full,  fitType);

% plot quality of fit.
hQualNocon = mctPlotQualityOfFit(r2_full,r2_orig, pSig_full,pSig_orig, dSig);

% Select the fibers explaining most of the varicance
hits_fibers = fiber_w >= 0.00001; % the fibers with large weights
fa_fibers  = ~hits_fibers;

fprintf('\n[%s] Extracting good and bad fibers...\n',mfilename)
fgGoodNocon = fgExtract(fgImg,find(hits_fibers),'keep');
fgBadNocon  = fgExtract(fgImg,find(fa_fibers),'keep');

% Plot fibers.
h1 = mctNfgDisplayStrands(fgImg,[.2 .75 .3],[],[],[],20,fiber_w(1:length(fgImg.fibers)) .* .1);
h2 = mctNfgDisplayStrands(fgGoodCon,[.75 .2 .3],[],[],[],20,.001.*ones(size(fgGoodCon.fibers)));
h3 = mctNfgDisplayStrands(fgBadCon,[.3 .2 .75],[],[],[],20,.001.*ones(size(fgBadCon.fibers)));

% Show the fibers
h4 = mctNfgDisplayStrands(fgLgnImg,[.2 .75 .3],[],[],[],20,fiber_w(1:length(fgLgnImg.fibers)) .* .1);
h5 = mctNfgDisplayStrands(fgGoodNocon,[.75 .2 .3],[],[],[],20,.001.*ones(size(fgGoodNocon.fibers)));
h6 = mctNfgDisplayStrands(fgBadNocon,[.3 .2 .75],[],[],[],20,.001.*ones(size(fgBadNocon.fibers)));

fileName  = sprintf('fe_test_lgn_or_%s_%s_algo%i',fitType,dtiDataType,algo);
saveMe = fullfile(saveDir,fileName);
fprintf('\n\nSaving file: %s\n\n',saveMe)
save(saveMe,'-v7.3')

saveas(h1,fullfile(saveDir,[fileName,'_fgAll_con']),'fig')
saveas(h2,fullfile(saveDir,[fileName,'_fgGood_con']),'fig')
saveas(h3,fullfile(saveDir,[fileName,'_fgBad_con']),'fig')

saveas(h4,fullfile(saveDir,[fileName,'_fgAll_nocon']),'fig')
saveas(h5,fullfile(saveDir,[fileName,'_fgGood_nocon']),'fig')
saveas(h6,fullfile(saveDir,[fileName,'_fgBad_nocon']),'fig')

saveas(hFitCon,fullfile(saveDir,[fileName,'_Fit_con']),'fig')
saveas(hQualCon,fullfile(saveDir,[fileName,'_Qual_con']),'fig')

saveas(hFitNocon,fullfile(saveDir,[fileName,'_Fit_nocon']),'fig')
saveas(hQualNocon,fullfile(saveDir,[fileName,'_Qual_nocon']),'fig')

keyboard

end


