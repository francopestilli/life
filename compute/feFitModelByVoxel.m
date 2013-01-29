function [fe voxfit w] = feFitModelByVoxel(fe)
%
% Fit the model one voxel at the time.
%
%  fefitvx = feFitModelByVoxel(fe)
%
% Franco (c) 2012 Stanford VISTA Team

%% Extract basic information regaridn gthe model.
nFibers = feGet(fe,'nfibers');
nVoxels = feGet(fe,'nvoxels');
 
tic
fprintf('[%s] Fitting the model voxel-wise...',mfilename)

%% Do the fit
% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

% Initialize the weights
w = nan(nFibers,nVoxels);
pSig = cell(nVoxels,1);
parfor ivox = 1:nVoxels
  [~, tmpw] = feFitModel(feGet(fe,'Mfiber',ivox), feGet(fe,'dsigdemeaned',ivox),'lsqnonneg');
  w(:,ivox) = tmpw';
  
  % Predict the signal.
  pSig{ivox} = feGet(fe,'Mfiber',ivox)*tmpw;
end

% Install the total weights and the predicted signal
w = w';
voxfit.weights = w;
voxfit.psig    = vertcat(pSig{:})';
fe             = feSet(fe,'voxfit',voxfit);

% Handling parallel processing.
if ~poolwasopen, matlabpool close; end
fprintf('done in %2.3fminutes.\n',toc/60)
      
return
