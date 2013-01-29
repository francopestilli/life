function s_mct_simulated_nfg
% function s_mct_simultated_nfg
%
% This function tests the MicroTrack related functions.
%
% More specifically, it tests the reduced and the full model on a small
% simulated dataset that.
%
% It simply shows that the fit of the reduced and full model can be represented as two separated linear models.
% The reduced model is fit to the de-meaned signal and generates correct predictions of the fiber weights.
% After this, the full model is fit to the signal (not de-meaned) and the
% estimated values for the baseline (isotropic) components are estimated
% correctly.
% 
% There are more tests we need to write and do beyound this. But this is a
% first solid results that shows that separating the full model into two
% components (fiber-weights and isotropic-weights) is meaningful. The full
% process of estimating the weights independently for fibers and weigths
% to then putting them back together to recover the original diffusion
% signal is not shown here but in test_mictrotrack_simulated_signal_recover.
%
% Example: 
%
% s_mct_simulated
%
% See also, test_mictrotrack_simulated_signal_recover
%
% (C) 2011 Stanford VISTA team. 


%% create a new fiber group uing nfg
% here we load a precomputed one
[fg dwi dwiNiftiFile] = mctNfgSimulate([],[],'strands_20120220T124436',[],[]);
nFiber         = length(fg.fibers);

nbvecs   =  size(find(dwi.bvecs(:,1)~=0),1);
bvecs    = dwi.bvecs; 
bvals    = dwi.bvals;
 

%% Make the list of tensors for each fiber and node
% These parameters could be adjusted to improve the quality of the fit.
d_ad = 1.5; d_rd = 0.5;
dParms(1) = d_ad; dParms(2) = d_rd; dParms(3) = d_rd;
fg.Q = fgTensors(fg,dParms);


%% Get the coordinates of the simulated fibers.
coords = fefgGet(fg,'unique image coords'); 


%% build the MicroTrack model
[Afiber, Aiso, dSig_simulated, dSig_simulated_demeaned, iso] = mctBuildDiffusionModel(dwi,fg,coords,[],'ones');
fprintf('\nThe A-redux matrix is %ix%i in size.\n',size(Afiber,1),size(Afiber,2))

%Afull = [Afiber Aiso];
%fprintf('\nThe A-full matrix is %ix%i in size.\n',size(Afull,1),size(Afull,2))

%% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data
fitType = 'lsqnonneg';

% This is a test for the L2 solution. If not many weights are 0 the L2
% solution should be correct within the noise added.
% fiber_w_l2 = Afiber \ dSig_simulated;
% assertAlmostEqual(fiber_w,fiber_w_l2,10^-2)

[fiber_w fiber_t]         = mctFitDiffusionModel(Afiber, dSig_simulated_demeaned, fitType);
fiber_s_w                 = mctSortModelWeights(fiber_w, nFiber);
[fiber_predDsig fiber_r2] = mctComputePredictedSignal(Afiber, fiber_w, dSig_simulated_demeaned) ;
mctDisplayModelFit(fiber_s_w, dSig_simulated_demeaned, fiber_t, fiber_predDsig, fiber_r2, fitType);


[full_w full_t]        = mctFitDiffusionModel(Afull, dSig_simulated, fitType);
full_s_w               = mctSortModelWeights(full_w, nFiber);
[full_predSig full_r2] = mctComputePredictedSignal(Afull, full_w, dSig_full_simulated) ;
mctDisplayModelFit(full_s_w, dSig_full_simulated, full_t, full_predSig, full_r2, fitType);

% Now let compute the weights for the isotropic components.
%
% These weights were implicitly set to the mean signal value in a voxel. 
% now that we solved for the fibers' weights we can estimate the isotropic
% weights.

% this is the residual signal we would liek to explain with the isotropic
% component. If we can think of the fibers explainign the variance of the
% de-meaned signal. This residual signal could be explained as the mean
% signal. 
residual_signal = dSig_simulated - fiber_predDsig;

% Build the Model (matrix) for the isotropic voxels only, this is a
% block-diagonal matrix.
Aiso = Afull(:,nFiber+1:end);
idx  = find(Aiso);
Aiso(idx) = 1;

% find the 'simple' isotropic component, the one we assumend for
% conscturcting the Afiber matrix, that is, the mean diffusion signal.
%
% We use an L2 minimization (least-square) to get the baseline signal. 
iso_w_L2 = Aiso \ residual_signal;

% Plot the results:
displayIsotropicFit(iso_w_L2, dSig_full_simulated, residual_signal, fiber_predDsig, nbvecs,coords)

end

% -------------------------------------------- %
function displayIsotropicFit(iso_w_L2, dSig_full_simulated, residual_signal, fiber_predDsig, nbvecs,coords)
   
mrvNewGraphWin('Isotropic component estimation'), 
% plot simulated signal and residual signal used to estimate the isotropic
% component (the mean of the signal in each voxel)
plot(residual_signal), hold on, 
plot(dSig_full_simulated,'r-');
plot(fiber_predDsig,'g-');
hold on
% plot he isotropic components on top of the signal used to estimate the
% to show that they are clsoe to the mean for the voxel they refer to:
x = nbvecs/2:nbvecs:size(coords,1)*nbvecs;
plot(x,iso_w_L2,'go','MarkerFaceColor','g','MarkerSize',15)
ylabel('DW signal')
legend( ...
sprintf('Residual signal used to estimate\nthe isotropic'), ...
        'Simulated signal', ...
        'Predicted isotropic signal', 'Location','NorthOutside') 

% mark each voxel on the plot
y = get(gca,'yLim');
y = repmat(y',1,size(coords,1));
x1 = x + x(1);
x1 = repmat(x1,size(y,1),1);
plot(x1, y,'k--')
for ii = 1:size(coords,1)
   voxLabel{ii} = sprintf('Voxel %i',ii); 
end

set(gca,'xTick',x(1,:), 'xTickLabel',voxLabel)
fprintf('The estimated isotropic components are: %0.2f\n',iso_w_L2)

end
