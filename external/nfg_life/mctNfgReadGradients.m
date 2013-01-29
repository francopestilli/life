function [bvecs bvals grad] = mctNfgReadGradients(nfgGradFile)
% 
% function fid = mctNfgReadGradients(nfgGradFile)
%
% Reads an NFG gradient file into bvecs and bvals
%
% Franco
%

grad = load(nfgGradFile,'-ascii');
bvecs = grad(:,1:3);
bvals = grad(:,4);


