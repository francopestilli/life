function fe = feRunLifeFit(fe)
% OBSOLETE
%
% fe = feRunLifeFit(fe,dwi)
%
% Performs a run of the life model and evaluates it against the
% original model
%
% Franco (C) 2012 Stanford VISTA team.

error('Obsolete %s\n',mfilename)

if ~isfield(fe,'life'), error('LiFE - missing ''life'' field in struct ''fe''.\n Run feBuildDiffusionModel.m first.');end

fprintf('\n[%s] the ROI used has size(%i,%i)\n',mfilename,size(fe.roi.coords))

% Fit the model.
fe = feFitConnectomeModel(fe);

% Fit the isotropic weights
% fe = feSet(fe,'isoweights',feGet(fe,'Aiso') \ feGet(fe,'dsig measured'));
% fe = feSet(fe,'fullweights',[feGet(fe,'fiber weights'); feGet(fe,'iso weights')]);

% Predict the signal
% fe = feSet(fe,'psig full',[]);
% fe = feSet(fe,'psig fiber',[]);

% Overall quality
% fe = feSet(fe,'r2 life',[]);
% fe = feSet(fe,'rmse',[]);

% Compute residuals, this is done on the demeaned (fiber) signal. 
% This signal was not explained by the model 
% fe = feSet(fe,'ressig',[]);

% Reshape the signals by voxel
% fe = feSet(fe,'psig vox',[]);
% fe = feSet(fe,'dsig vox',[]);
% fe = feSet(fe,'ressig vox',[]);

% Compute R2 and rmse by voxel.
% fe = feSet(fe,'r2 vox',[]);
% fe = feSet(fe,'rmse vox',[]);

% Fiber density, fiber count etc??
% FIX FIX

end