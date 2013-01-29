function [fascicles,classification,fg,fg_classified] = feAfqSegment(dtFile, wholeBrainConnectome)
% This runs an AFQ segmenation on a whole brain fiber group. Clean the
% fiber groups and return indices to the fascicles can can be used to
% address the fibers inside LiFE.
% 
%    [fascicles,classification,fg,fg_classified] = feAfqSegment(dtFile, wholeBrainConnectome)
%
% INPUTS:
%  dtFile               -  fullpath to a dt6 file 
%  wholeBrainConnectome - the full path to a whole-brain connectome  
%
% OUTPUTS:
%  fascicles      - tCell array with the 20 fascicles classified by AFQ
%  fg             - the whole brain conenctome cleaned by AFQ, these are all the fibers
%                   use by AFQ to create the 20 fascicles. THese fibers are less than the
%                   in the original connectome.
%  classification - a vector of indices indicating to which fascicle in fg_classified each
%                   fiber in fg is part of.
%  AFQ 'classified_fg' - a single fiber group with all the fascicles
%                        indexed by the subfiled subgroup.
%
% Franco (C) 2012 Stanford VISTA team.

% Load the connectome
fg = fgRead(wholeBrainConnectome);

% Segment the fibers using AFQ
[fg_classified,~,classification,fg]= AFQ_SegmentFiberGroups(dtFile, fg);

% Split the fiber groups into individual groups
fascicles = fg2Array(fg_classified);

return