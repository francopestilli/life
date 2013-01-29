function params = mctNfgReadParams(nfgParamFile)
% 
% function fid = mctNfgReadParams(nfgParamFile)
%
% Reads an NFG parameters file into a structure.
% params{1} shows the parameter ID
% params{2} shows the parameter VALUE
%
% Values are loaded as strings, because some parameters are strings and
% others are integers.
%
% Franco
%
fid = fopen(nfgGradientFile,'rt');
params = textscan(fid,'%s%s','space');
fclose(fid);