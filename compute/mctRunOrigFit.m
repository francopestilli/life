function fe = mctRunOrigFit(fe,dwi)
%
% fe = mctRunOrigFit(fe,dwi)
%
% Performs a run of the orginal connectome model.
%
% The original model is a model that assigns the same weight to all the
% fibers in the connectome.
%
% Franco (c) 2012 Stanford VISTA team.

if ~isfield(fe,'life')
   error('This function requires to run first: fe = mctRunLifeFit(fe);')
end

% Compute the diffusion signal predicted by the original fiber model.
AorigFull  = [sum(fe.life.Afiber,2)  fe.life.Aiso * fe.life.w.iso]; 
w          = AorigFull \ fe.dSig;

fe.orig.w.fiber  = w(1);
fe.orig.w.iso    = w(2);
fe.orig.w.full   = w;
fe.orig.pSigFull = mctComputePredictedSignal(AorigFull,fe.orig.w.full);
fe.orig.pSig     = mctComputePredictedSignal(AorigFull(:,1),fe.orig.w.fiber);

% Overall quality across voxels
[fe.orig.rmse, fe.orig.r2] = mctComputePredictionQuality(fe.dSig_demeaned, fe.orig.pSig,2);

% Compute residual signal
fe.orig.resSig = (fe.dSig - fe.orig.pSig) + AorigFull(:,2);

% Reshape the signal by voxel and compute qulity of fit by voxels
fe.orig.vox.pSig   = mctReshape(full(fe.orig.pSig), fe.nBvecs, fe.nVoxels);
fe.orig.vox.pSigUW = mctReshape(full(fe.orig.pSig) * (1/fe.orig.w.fiber), fe.nBvecs, fe.nVoxels);
fe.orig.vox.dSig   = mctReshape(full(fe.dSig_demeaned), fe.nBvecs, fe.nVoxels);
[rmse, r2]         = mctComputePredictionQuality(fe.orig.vox.dSig, fe.orig.vox.pSig,2);
fe.orig.vox.r2     = r2;
fe.orig.vox.rmse   = rmse;
fe.orig.vox.res    = mctReshape(full(fe.orig.resSig), fe.nBvecs, fe.nVoxels);
