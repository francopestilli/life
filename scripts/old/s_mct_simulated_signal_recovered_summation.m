function test_microtrack_simulated_signal_recovered_summation
% function  test_microtrack_simulated_signal_recovered_summation
%
% This function tests the MicroTrack related functions.
%
% More specifically, it tests that when the signal is simulated given two
% fibers but attempted to be recovered with onlyone of the two fibers the
% recovery is suboptimal.
%
% Example: test_microtrack_simulated_signal_recovered_summation
%
% See also: test_microtrack_simulated_recovered
%
% (C) 2011 Stanford VISTA team.
%
% To do: Make example in which one fascicle is removed. Out of the data and
% the whole data set is attemtped to be explained by both.

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
fg.pathwayInfo = {{},{}};

%% build the MicroTrack model
[Afull dSig] = mctBuildFullDiffusionModel(dwi,fg,coords);
fprintf('\nThe A-full matrix is %ix%i in size.\n',size(Afull,1),size(Afull,2))


%% now set some weights and compute dSig with those weights
noiseLevel = 0.2; % this is the noise level relative to the signal

% full model created with two fibers:
sim_weights = [.9 .1 ones(1,size(coords,1))]';
sim_dSig    = Afull * sim_weights;

% corrupt with noise
noiseSd = noiseLevel * std(sim_dSig);
noise   = noiseSd   .* randn(size(sim_dSig));

% this si the signal that we simulate and that we will want to use to
% estimate the fiber and the isotropic weights BUT also this is the signal
% that will want to recover at the end of the whole busisness.
%
% We will attempt to recover this signal by combinign the weights estimated
% on the fibers and on the isotropic components independently.
sim_dSig_noise = noise + sim_dSig;


%% Find the fibers weights.

% chose a a fitting engine.
fitType = 'cvx';

% Find the smallest number fo fibers that explain most of the variance in
% the diffusion data.

% Note, the matrix here is created with only one fiber, but the signal was
% created with two fibers.
fg1 = fgExtract(fg,1);

keyboard

% this currently does nto make sense.
% because removing one fiber from the group implies removing the unique
% voxels that fiber goes through, which in turns implies changing the A
% matrix. Which properly hat we do not want.
% THINK MORE.

% we have to subselect the voxels to use for the model
coords1 = fefgGet(fg1,'unique image coords'); 
[Afiber dSig, dSig_demeaned iso] = mctBuildDiffusionModel(dwi,fg1,coords1);

% fit the model with only one fiber on a signal created with two fibers.
[fiber_w fiber_t]         = mctFitDiffusionModel(Afiber, sim_dSig_noise, fitType);

% show the results of the fit. Please note that the scatter plot is not
% expected to look good, this is because we used a de-meaned Afiber model
% with a signal with the mean value in it. So we cannot really account for
% the DC components.
w_sorted            = mctSortModelWeights(fiber_w, nFiber);

% predict the signal from the fiber weights. Note, the above note abotu signal mean :) 
[fiber_predDsig fiber_r2] = mctComputePredictedSignal(Afiber, fiber_w, sim_dSig_noise) ;
mctDisplayModelFit(w_sorted, sim_dSig_noise, fiber_t, fiber_predDsig, fiber_r2, fitType);

%% Find the weights for the isotropic (DC) components

% These weights were implicitly set to the mean signal value in a voxel. 
% now that we solved for the fibers' weights we can estimate the isotropic
% weights.

% We first compute the residual signal. This is the signal we would like to
% explain with the isotropic component. That is the signal that is not
% explained by the fiber weights. If we can think of the fibers explainign
% the variance of the de-meaned signal. This residual signal could be
% explained as the mean signal.
residual_signal = sim_dSig_noise - fiber_predDsig;

% Now we build the model for the isotropic components.
% This is the block-diagonal matrix built for the right-hand sida of Afull.
Aiso = Afull(:,nFiber+1:end);
idx  = find(Aiso);
Aiso(idx) = 1; % we put ones here for now because we are not tinking about units for the moment

% Find the weights for isotropic (DC) components.
% We use an L2 minimization (least-square) to get this baseline signal. 
tic;
iso_w = Aiso \ residual_signal;
% iso_w = iso; % let think about why this line can be substituted for the
% precedent and everything works
t = toc;

% Plot the results:
displayIsotropicFit(iso_w, sim_dSig_noise, residual_signal, fiber_predDsig, nbvecs,coords)

%% Now build the full model, using the fiber and isotropic weights estimated independently.

% model weights, fibers + isotropic components (arbitrary units).
w = [fiber_w', iso_w']';

% this is the final model, where the fibers estimators are de-meaned but
% the isotropic components estmators are not (they are set to 1's).
A = [Afiber, Aiso];

% Replot the model fit after having estimated the weights for the isotropic
% component.
[recovered_dSig recovered_r2] = mctComputePredictedSignal(A, w, sim_dSig_noise) ;
w_sorted.isotropic = iso_w;
mctDisplayModelFit(w_sorted, sim_dSig_noise, fiber_t + t, recovered_dSig, recovered_r2, fitType);



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
        'Fuber predicted signal (de-meaned)', ...
        'Predicted isotropic signal (DC)', 'Location','NorthOutside') 

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
