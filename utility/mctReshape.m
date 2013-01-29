function val = mctReshape(val, nBvecs, nVoxels)
%
% function val = mctReshape(val, nBvecs, nVoxels)
%
% THis function reshapes the signals from vector form to image (volume)
% form.
%
% This simple operations is called several times to transform the vectors
% of signal, residuals etc into t a shape that can be saved as an ROI and
% displayed on the tstructure in mrDiffusion.
%
%
% See also:  mctComputePredictionQuality.m
%
% Example:
%  val = mctReshape(randn(150*57,1), 150, 57)
% 
% (c) Stanford VISTA Team

val = reshape(val,nBvecs,nVoxels);
