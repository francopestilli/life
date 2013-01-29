function s_mct
% function s_mct
% Some microtrack analysis (arcuate fasciculus). 
% 
% (C) 2011 Stanford VISTA team. 


%% Read the fibers and dwi data.
dataDir = fullfile(mrvDataRootPath,'diffusion','fiberPrediction');
dwi     = dwiLoad(fullfile(dataDir,'raw','dti_g13_b800_aligned.nii.gz'));

% Fibers are stored in ACPC space.
% We load up the dt6 and open a mrDiffusion window.  This brings in the
% xforms so that we can easily transform between image and acpc space.
dt6Name      = fullfile(dataDir,'dti06','dt6.mat');
[dtiF, dtiH] = mrDiffusion('off',dt6Name);

%% XXX Need to work here on the spatial transformations: 
fgName = fullfile(mrvDataRootPath,'diffusion','sampleData','fibers','leftArcuate.pdb');

% get te xForms, from a dn to Acpc
xForm.img2acpc = dtiGet(dtiH,'img 2 acpc xform');
xForm.acpc2img = dtiGet(dtiH,'acpc 2 img xform');

fgAcpc = mtrImportFibers(fgName,xForm,[],'acpc'); % Fiber coordinates in acpc space
fgImg  = dtiXformFiberCoords(fgAcpc,xForm.acpc2img,'img');

%% extract a subset of the fibers so to have a better handle on them.
fibersToKeep = 232; % number of fibers to keep and use for fitting and testing.

% we keep small fibers, with less than 80 nodes, this will make A small.
fibersNodes  = find(fefgGet(fgImg,'nodes per fiber') < 100); 
fg           = fgExtract(fgImg, randsample(fibersNodes,fibersToKeep),'keep');
nFiber       = fefgGet(fg,'n fibers');

%% get the coordinates the fibers go through 
coords = fefgGet(fg,'unique image coords'); 

%% Make the list of tensors for each fiber and node
% These parameters could be adjusted to improve the quality of the fit.
d_ad = 1.5; d_rd = 0.5;
dParms(1) = d_ad; dParms(2) = d_rd; dParms(3) = d_rd;
fg.Q = fgTensors(fg,dParms);

%% build the MicroTrack model
% Afiber is the portion for the fibers
% Aiso is the isotropic portion
% dSig is the diffusion signal (raw)
% dSig_demeaned is the diffusion signal with the mean removed.
[Afiber, Aiso, dSig, dSig_demeaned] = mctBuildDiffusionModel(dwi,fg,coords,[],'ones');

%% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data

% cvxPath - on Windows for Wandell
fitType = 'tfocsl2';  

[fiber_w fiber_t]   = mctFitDiffusionModel(Afiber, dSig, fitType);
keyboard

w_sorted            = mctSortModelWeights(fiber_w, nFiber);
[fiber_predDsig fiber_r2] = mctComputePredictedSignal(Afiber, fiber_w, dSig);

% plotting fiber predictions against the signal  
mctDisplayModelFit(w_sorted, dSig, fiber_t, fiber_predDsig, fiber_r2, fitType);

% plotting the fiber prediction against the de-meaned signal
mctDisplayModelFit(w_sorted, dSig_demeaned, fiber_t, fiber_predDsig, fiber_r2, fitType);

residual_signal = dSig - fiber_predDsig;


% Now we build the model for the isotropic components.
% This is the block-diagonal matrix built for the right-hand sida of Afull.

% Find the weights for isotropic (DC) components.
% We use an L2 minimization (least-square) to get this baseline signal. 
tic;
iso_w = Aiso \ residual_signal;
% iso_w = iso; % let think about why this line can be substituted for the
% precedent and everything works
t = toc;

%% Now build the full model, using the fiber and isotropic weights estimated independently.

% model weights, fibers + isotropic components (arbitrary units).
w = [fiber_w', iso_w']';

% this is the final model, where the fibers estimators are de-meaned but
% the isotropic components estmators are not (they are set to 1's).
A = [Afiber, Aiso];

% size(A)

% Replot the model fit after having estimated the weights for the isotropic
% component.
[recovered_dSig recovered_r2] = mctComputePredictedSignal(A, w, dSig) ;
w_sorted.isotropic = iso_w;
mctDisplayModelFit(w_sorted, dSig, fiber_t + t, recovered_dSig, recovered_r2, fitType);

% Show the fibers
% [~, goodFibers] = find(fiber_w > 0.0251);
% fgView(dtiH,fgExtract(fg,goodFibers),175);

keyboard

end

function matrix_size(A)

fprintf('\nThe A matrix is %ix%i in size.\n',size(A,1),size(A,2))

end