function s_life_rois_wholebrain(fitType, lambda,algo, this_roi)
% function s_life_rois_wholebrain
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
if notDefined('this_roi');   this_roi = 'mct_r_lgn.mat';end
if notDefined('fitType');   fitType = 'sgd';end

dataRootPath = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
subfolders   = fullfile('150dirs_b1000_1');
dataDir      = fullfile(dataRootPath,subfolders);
dtFile       = fullfile(dataDir,'dt6.mat');
dwiFile      = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');

% PRecomputed ROIs and Fiber groups
roiDir       = fullfile('/biac2/wandell6/data/arokem/ModelFits/FP20120420/',subfolders,);
fiberDir     = fullfile('/biac2/wandell6/data/arokem/ModelFits/FP20120420/',subfolders,'fibers','whole_brain_FACT');
fgName       = fullfile(fiberDir,'WholeBrainFG_FACT.mat');

% load the dwi and dti files        
dwi          = dwiLoad(dwiFile);
dt           = dtiLoadDt6(dtFile);
[dtiF, dtiH] = mrDiffusion('off',dtFile);
 
% get te xForms, from a dn to Acpc
xForm.img2acpc = dtiGet(dtiH,'img 2 acpc xform');
xForm.acpc2img = dtiGet(dtiH,'acpc 2 img xform');
close(dtiF); clear dtiH dtiF


%% Get the ROI
roiName    = fullfile(roiDir,this_roi);
roi        = dtiReadRoi(roiName);
roi.coords = floor(mrAnatXformCoords(xForm.acpc2img ,roi.coords));
coords     = unique(roi.coords,'rows');
fprintf('\n[%s] ROI: %s, size(%i,%i)\n',mfilename,roi.name,size(roi.coords))


%% Load the WHOLE BRAIN fibers
fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,fgContrack)
fg     = dtiLoadFiberGroup(fgName);
fgImg  = dtiXformFiberCoords(fg, xForm.acpc2img,'img');
clear fg

% Select only the fibers' coordinates going throguht the ROI
% and only the ROI coordinates where a fiber goes through
unique_fg_coords = fefgGet(fgImg,'unique image coords');
coords = intersect(coords,unique_fg_coords,'rows');
fprintf('\n[%s] ROI IMG %s, size(%i,%i)\n',mfilename,roi.name,size(coords))

% Clip the fibers to the portion inside the ROI
fg = dtiIntersectFibersWithRoi([],[],[],roi,fgImg)
keyboard

%% Build the MicroTrack model
% Afiber is the portion for the fibers
% Aiso is the isotropic portion
% dSig is the diffusion signal (raw)
% dSig_demeaned is the diffusion signal with the mean removed.
[Afiber, Aiso, dSig, dSig_demeaned, dSig_fiberPredicted] = mctBuildDiffusionModel(dwi,fgImg,coords,[],'ones');
nFiber = size(Afiber,2); 

% Check that none of the columns has all-zeros, meaning some fibers do not
% go throguht this ROI.
% all(Afiber==0,1)

% build the full model
Afull = [Afiber, Aiso];

%% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data
% Fit the fibers
fiber_w = mctFitDiffusionModel(Afiber, dSig_demeaned, fitType,lambda);

% Fit the isotropic
iso_w = Aiso \ dSig;
full_w = [fiber_w'; iso_w];

% save results to disk
saveDir = '/biac2/wandell6/data/frk/LiFE';
fileName  = sprintf('%s_%s_%s_lambda%i_algo%i',this_roi(1:end-4),fitType,dtiDataType,10000*lambda,algo);
saveMe = fullfile(saveDir,fileName);
fprintf('\n\nSaving file: %s\n\n',saveMe)
save(saveMe,'-v7.3')

%% Predict the signal 
pSig_full  = mctComputePredictedSignal(Afull,full_w);
pSig_fiber = mctComputePredictedSignal(Afiber,fiber_w);

%% Compute prediction quality of LiFE WM model
nVoxels = size(Aiso,2); 
nBvecs  = size(Afiber,1) / nVoxels;

% Overall quality across voxels
[rmse_full, r2_full] = mctComputePredictionQuality(dSig, pSig_full);

% quality by voxel
%dSig_img      = mctReshape(dSig, nBvecs, nVoxels);
%pSig_full_img = mctReshape(pSig_full, nBvecs, nVoxels);
%[rmseVox_full, r2vox_full] = mctComputePredictionQuality(dSig_img, pSig_full_img);

% Overall quality across voxels
[rmse_fiber, r2_fiber] = mctComputePredictionQuality(dSig_demeaned,pSig_fiber);

% Quality by voxel
%dSig_dem_img      = mctReshape(dSig_demeaned, nBvecs, nVoxels);
%dSig_fiber_img = mctReshape(pSig_fiber, nBvecs, nVoxels);
%[rmseVox_fiber, r2vox_fiber] = mctComputePredictionQuality(dSig_dem_img, dSig_fiber_img);

% Compute the fibers density in each voxel. 
%[fiberDensity_life, nodesNum_life] = mctComputeFiberDensity(fgImg,coords,fiber_w);
%ratioOfNodesNum_life = (max(nodesNum_life) / min(nodesNum_life));


%% Compute prediction quality of the original WM model
% Compute the diffusion signal predicted by the original fiber model.
Aorig     = [sum(Afiber,2)  Aiso * iso_w]; 
orig_w    = Aorig \ dSig;
pSig_orig = mctComputePredictedSignal(Aorig,orig_w);

% Overall quality across voxels
[rmse_orig, r2_orig] = mctComputePredictionQuality(dSig, pSig_orig);

% quality by voxel
%pSig_orig_img = mctReshape(pSig_orig, nBvecs, nVoxels);
%[rmseVox_orig, r2vox_orig] = PmctComputePredictionQuality(dSig_img, pSig_orig_img);

% Compute the fibers density in each voxel. 
%[fiberDensity_orig, nodesNum_orig] = mctComputeFiberDensity(fgImg,coords);
%ratioOfNodesNum_orig = (max(nodesNum_orig) / min(nodesNum_orig));


%% Plot the results
% parse fiber weights and isotropic weights
[wFibers_fu wIsotropic] = mctSortModelWeights(full_w,nFiber);
wFibers_fi = mctSortModelWeights(fiber_w,nFiber);

% plotting fiber predictions against the signal  
h1 = mctDisplayModelFit(wFibers_fu,wIsotropic, dSig, pSig_full, r2_full,  fitType);

% plot quality of fit.
h2 = mctPlotQualityOfFit(r2_full,r2_orig, pSig_full,pSig_orig, dSig);


%% Select the fibers explaining most of the varicance
hits_fibers = fiber_w >= 0.00001; % the fibers with large weights
fa_fibers  = ~hits_fibers;

fgGood = fgExtract(fgImg,find(hits_fibers),'keep');
fgBad  = fgExtract(fgImg,find(fa_fibers),'keep');

% save results to disk
fprintf('\n\nSaving file: %s\n\n',saveMe)
save(saveMe,'-v7.3')
savefigvista(h1,[fileName,'_weigths'],'eps',saveDir,'/',1,1);
savefigvista(h2,[fileName,'_r2'],'eps',saveDir,'/',1,1);

%% Write down the results
%mctWriteStats2nifti(fiberDensityAfter,fullfile(dataDir,'fiber_density_after'),dt,coords);
%mctWriteStats2nifti(fiberDensityBefore,fullfile(dataDir,'fiber_density_before'),dt,coords);

%% Show the seelcted fibers
%fgView(dtiH,fgGood,300);
%fgView(dtiH,fgBad,301);
%fgView(dtiH,fg,302)

end

%% test lambda values
% 
% lambda = [.001 .02 1 20 100];
% for ii = 1:length(lambda)
%   close all;
%   clc;
%   s_life_rois_wholebrain(lambda(ii),'mct_prefrontal_10mm.mat');
% end