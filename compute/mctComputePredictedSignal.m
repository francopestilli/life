function predSig = mctComputePredictedSignal(A,w) 
%
% function predSig = mctComputePredictedSignal(A,w) 
%
% Compute the predicted signal from the microtrack model
%
%
%  See also s_mct_fact_cc_prefrontal_rois.m
%
%  Example:
%    predSig = mctComputePredictedSignal(A,w) 
%
% (c) Stanford VISTA Team

% Transponse the weights if they were passed with the wrong dimensions
if ~(size(A,2) == size(w,1)), w =w';end
predSig = A * w;

return

