function s_mct_simulated
% function s_mct_simultated
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


%% create a new fiber group
crossing = 1;
fg = fgSimulateHelper([],crossing);
nFiber         = length(fg.fibers);


%% create a new dwi structure
brainSiz = 100;
nbvecs   = 100; 
meanDiffusivity = 1000;
noiseSd         = 0.2;

% set data to a mean diffusivity
data     = meanDiffusivity + zeros(brainSiz, brainSiz,brainSiz,  nbvecs + 1);

% add noise to data
data = data + (noiseSd .* randn(size(data)));
nifti    = niftiGetStruct(data);
nifti    = setfield(nifti,'descrip','Simulated diffusion nifti data');

% Make bvecs with b0 on top of it: 
bvecs    = [makeBvecsHelper(nbvecs); 0 0 0]; 

% bvals are generally refered to in the hundreds-thousands but here we have them coded in as smaller, why? 
bValue = 900;
bvals  = [bValue/1000 * ones(nbvecs,1); 0];
dwi    = dwiCreate('name','Simulated DWI data','nifti',nifti,'bvecs',bvecs,'bvals',bvals);
 

%% Get the coordinates of the simulated fibers.
coords = fefgGet(fg,'unique image coords'); 


%% Make the list of tensors for each fiber and node
% These parameters could be adjusted to improve the quality of the fit.
d_ad = 1.5; d_rd            = 0.5;
dParms(1) = d_ad; dParms(2) = d_rd; dParms(3) = d_rd;
fg.Q = fgTensors(fg,dParms);


%% build the MicroTrack model
[Aredux dSig, dSig_demeaned iso] = mctBuildDiffusionModel(dwi,fg,coords);
fprintf('\nThe A-redux matrix is %ix%i in size.\n',size(Aredux,1),size(Aredux,2))

[Afull dSig] = mctBuildFullDiffusionModel(dwi,fg,coords);
fprintf('\nThe A-full matrix is %ix%i in size.\n',size(Afull,1),size(Afull,2))


%% now set some weights and compute dSig with those weights
noiseLevel = 0.5; % this is the noise level relative to the signal

% reduced model (fiber prediction only)
weights_redux        = [.9 .5]';
dSig_redux_simulated = Aredux * weights_redux;

% corrupt with noise
noiseSd = noiseLevel * std(dSig_redux_simulated);
noise   = noiseSd   .* randn(size(dSig));
dSig_redux_simulated_noise = noise + dSig_redux_simulated;

% full model (fibers + isotropic component)
weights_full        = [.9 .5 ones(1,size(coords,1))]';
dSig_full_simulated = Afull * weights_full;

% corrupt with noise
noiseSd = noiseLevel * std(dSig_full_simulated);
noise   = noiseSd   .* randn(size(dSig));
dSig_full_simulated_noise = noise + dSig_full_simulated;

% here we show that the signals were corrupted by noise to a similar extent:
fprintf('Correlation between noise- and not noise-corrupted signal\nfor the FULL  model r2=%0.3f and for the REDUX model r2=%0.3f\n', ...
         corr2(dSig_full_simulated_noise,dSig_full_simulated)^2, ...
         corr2(dSig_redux_simulated_noise,dSig_redux_simulated)^2)
disp('If these two r2''s are very similar signals were corrupted with similar relative noise.\n\n')


%% Fit the model.
% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data
fitType = 'lsqnonneg';

% This is a test for the L2 solution. If not many weights are 0 the L2
% solution should be correct within the noise added.
% redux_w_l2 = Aredux \ dSig_redux_simulated;
% assertAlmostEqual(redux_w,redux_w_l2,10^-2)

[redux_w redux_t]         = mctFitDiffusionModel(Aredux, dSig_redux_simulated_noise, fitType);
redux_s_w                 = mctSortModelWeights(redux_w, nFiber);
[redux_predDsig redux_r2] = mctComputePredictedSignal(Aredux, redux_w, dSig_redux_simulated_noise) ;
mctDisplayModelFit(redux_s_w, dSig_redux_simulated_noise, redux_t, redux_predDsig, redux_r2, fitType);

% This is a test for the L2 solution. If not many weights are 0 the L2
% solution should be correct within the noise added.
% full_w_l2 =  Afull \ dSig_full_simulated;
% assertAlmostEqual(full_w,full_w_l2,10^-2)

[full_w full_t]        = mctFitDiffusionModel(Afull, dSig_full_simulated_noise, fitType);
full_s_w               = mctSortModelWeights(full_w, nFiber);
[full_predSig full_r2] = mctComputePredictedSignal(Afull, full_w, dSig_full_simulated_noise) ;
mctDisplayModelFit(full_s_w, dSig_full_simulated_noise, full_t, full_predSig, full_r2, fitType);

% Now let compute the weights for the isotropic components.
%
% These weights were implicitly set to the mean signal value in a voxel. 
% now that we solved for the fibers' weights we can estimate the isotropic
% weights.

% this is the residual signal we would liek to explain with the isotropic
% component. If we can think of the fibers explainign the variance of the
% de-meaned signal. This residual signal could be explained as the mean
% signal. 
residual_signal = dSig_full_simulated_noise - redux_predDsig;

% Build the Model (matrix) for the isotropic voxels only, this is a
% block-diagonal matrix.
Aiso = Afull(:,nFiber+1:end);
idx  = find(Aiso);
Aiso(idx) = 1;

% find the 'simple' isotropic component, the one we assumend for
% conscturcting the Aredux matrix, that is, the mean diffusion signal.
%
% We use an L2 minimization (least-square) to get the baseline signal. 
iso_w_L2 = Aiso \ residual_signal;

% Plot the results:
displayIsotropicFit(iso_w_L2, dSig_full_simulated_noise, residual_signal, redux_predDsig, nbvecs,coords)

end

% -------------------------------------------- %
function displayIsotropicFit(iso_w_L2, dSig_full_simulated_noise, residual_signal, redux_predDsig, nbvecs,coords)
   
mrvNewGraphWin('Isotropic component estimation'), 
% plot simulated signal and residual signal used to estimate the isotropic
% component (the mean of the signal in each voxel)
plot(residual_signal), hold on, 
plot(dSig_full_simulated_noise,'r-');
plot(redux_predDsig,'g-');
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

% fgSimulateHelper &
function fg = fgSimulateHelper(fg,crossing,varargin)
%
% This helper function gets' a fiber group and changes the fibers to a
% hard-coded simple example with either 2 paraller fibers (crossing=0) or 2
% crossing fibers (crossing=1)
%
% 
% franco

if isempty(fg)
    fg = fgCreate('name','Simulated fiber group','coordspace','img');
end

switch crossing    
    case {1} % two crossing fibers
        fg.name = 'Two simulated crossing fibers';
        fg.colorRgb  = [200 50 50];
        fg.fibers    = {};
        fg.fibers{1} = [2,2,1; 3,2,1; 4,2,1]';
        fg.fibers{2} = [3,1,1; 3,2,1; 3,3,1]';
        
    case {0} % two parallel fibers
        fg.name = 'Two simulated parallel/independent fibers';
        fg.colorRgb  = [50 200 50];
        fg.fibers    = {};
        fg.fibers{1} = [2,2,1; 3,2,1; 4,2,1; 5,2,1]';
        fg.fibers{2} = [2,3,1; 3,3,1; 4,3,1; 5,3,1]';
    
    otherwise
     % crossing is used a s a parameter the set the total number of fibers
     fg.name = sprintf('%i simulated fibers',crossing);
     fg.colorRgb  = [50 50 200];
     fg.fibers    = {};

     if ~isempty(varargin), brainSiz = varargin{1};end
     avrgFiberLength = ceil(0.3 * brainSiz); % fibers' length is 30% of brain
     
     for ii = 1:crossing
         % FIX FIX FIX randperm(brainSiz) how to do this!
         fg.fibers{ii} = [2,2,1; 3,2,1; 4,2,1; 5,2,1]';
     end
        
end
end


% makeBvecsHelper %
function bvecs = makeBvecsHelper(n)
% 
% Returns the cartesian coordinates of n bvecs which are maximally distant
% on the surface of the sphere, computed with an electrostatic repulsion
% model
% 

bd_dir = sprintf('%s/mrDiffusion/data/gradFiles/',mrvDirup(mrvRootPath,1)); 

pts_dir = fullfile(bd_dir,'caminoPts'); 

elec_points = dlmread(fullfile(pts_dir, sprintf('Elec%03d.txt', n))); 

% The first line should be equal to n: 
assert(elec_points(1)==n, 'There is something wrong with the camino points file'); 

% The format is: n,x1,y1,z1,x2,y2,z2 ... xn,yn,zn 
bvecs = [elec_points(2:3:end),... 
               elec_points(3:3:end),...    
               elec_points(4:3:end)];
end
