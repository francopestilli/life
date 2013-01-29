function [mad, r2] = mctComputeDataReliability(dSig)
% 
% function [] = mctComputeDataReliability(dSig)
%
% Compute the reliability of the data assuming that several repetitions of
% the same measurement is available.
%
% See also: s_mct_fact_cc_prefrontal_rois.m, mctComputeDataReliability,
%           mctDiffusionModel.m
%
% Example:
%  See [mad, r2] = mctComputeDataReliability(dSig)
% 
% (c) Stanford VISTA Team

% ideally i would like to use different cross validation methods.
% for the moment I implement a test-retest method.

sig_test   = dSig(:,1:floor(size(dSig,2)/2));
sig_retest  = dSig(:,floor(size(dSig,2)/2)+1:end);

[mad, r2] = mctComputePredictionQuality(sig_test,sig_retest);

return
