function [fg_classified,classification,fg] = s_life_afq(dtFile, wholeBrainConnectome)
% function s_life_afq
%
% This runs an AFQ segmenation on a whole brain fiber group.
%
% INPUT:
%  dwiData - the full path to a  
%
% Franco
%
% (C) 2012 Stanford VISTA team.

%% Load the connectome
fg = fgRead(wholeBrainConnectome);

%% Segment the fibers using AFQ
[fg_classified,~,classification,fg]= AFQ_SegmentFiberGroups(dtFile, fg);

return
