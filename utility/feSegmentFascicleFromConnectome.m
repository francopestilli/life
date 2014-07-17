function [fg, keepFG] = feSegmentFascicleFromConnectome(fg, rois, operation, fascicleFileName)
% Segment a fascicle from a connectome by applying a series of 'AND' and
% 'NOT' operations between the connectome and a set of ROIs.
%
% The important aspect of this function is that it returns both the new
% fiber group AND the indices in the origianl fiber group of the fibers
% that passed all the 'and' and 'not' operations requested.
%
% These idices can be used to address the columns of a LiFE model (M
% matrix). This allows for subtracting entire fiber group dfined
% anatomically from a pre-existing connectome. In turns this allows for
% testing hypotheses on the improtance of an (anatomically selected)
% fascicle within the volume of white-mater comprised by the connectome.
%
% [fg fibersIndices] = feSegmentFascicleFromConnectome(fg, rois, operation, fascicleFileName)
%
% INPUTS:
%   fg           - A connectome (e.g., fg = feGet(fe,'fibers acpc'))
%   rois         - A cell array of rois to be used for selecting the fibers in the
%                  final fascicle
%   operation    - A set of logical operations to be applied to the fibers in
%                 the connectome in relation to the rois. There should be one
%                 operation per roi. 
%   fascicleName - The name of the final fascicle. 
%
%
% OUTPUTS:
%  fg            - The segmented fiber group, containing only the fibers
%                  that passed all the logical operations.
%  keepFG        - A vector of indices (0's and 1's) of the length of the
%                  number of fibers in the input fiber group. A one
%                  indicates that the fiber survived all the logical
%                  operations requested. A 0 indicates that the fibers did
%                  not survive some of the logical operations and it was
%                  deleted from the output fiber group.
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com.

% Make sure that the inputs have one logical operation per ROI
if ~(length(rois) == length(operation))
  error('[%s] Please provide one logical operand (e.g., ''and'', ''not'') for each ROI...',mfilename)
end

% Read fibers, if a path was passed
if ~isstruct(fg)
  fg = fgRead(fg);
end

% The following is the vector containing the indices 
% to the fibers we KEEP from the origianl fiber group
% after all the logical operations are applied.
keepFG = false(length(fg.fibers),1); 

% The following is a cell array which will old the relative indices of the
% fibers into each sized-down version of the fibers.
% Each entry of the cell arry holds the indices to the fibergroup in the
% before the current operation was applied.
currentFibIndices = cell(length(rois)+1,1);    
currentFibIndices{1} = 1:length(keepFG);

for ir = 1:length(rois)
  % Read the rois from disk if paths were passed in
  if ~isstruct(rois{ir})
    rois{ir} = dtiReadRoi(rois{ir});
  end
  
  % Intersect the wholebrain fiber group with "AND" / "NOT" ROIs
  [fg, ~, keep]  = dtiIntersectFibersWithRoi([],operation{ir},[],rois{ir},fg);
  
  % Select the indices fo the fibers that were deleted in the previous
  % loop. The way we address these indices depends on the type of operation.
  % For 'and' we simply use the idices in keep{ir-1}
  % For 'not' we need to flip the sign and use the indices in kee{ir-1}
  % that were se to 0.
  %
  % See help dtiIntersectFibersWithRoi.m
  % "The output variables keep and keepID vectors for "not" option are
  % counterintuitive: they mark fibers that DO intersect the ROI and that are
  % exluded from the output FG."
  switch operation{ir}
    case {'and','AND','and both endpoints','endpoints'}
      currentFibIndices{ir+1} = currentFibIndices{ir}(keep);
    case {'not','NOT'}
      currentFibIndices{ir+1} = currentFibIndices{ir}(~keep);
    otherwise
      keyboard
  end
  clear keep
end

% Save the indices of the fibers that survived all the operations.
keepFG( currentFibIndices{end} ) = true;

% Clean up the rest of the fields in the fiber group.
% dtiIntersectFibersWithRoi.m does not handle other fields but the .fiber
% one.
fg.pathwayInfo = [];
fg.seeds       = [];
fg.Q           = [];
fg.params      = [];

% Change the fibergroup name
if isstruct(fascicleFileName)
    fg.name = fascicleFileName.name;
else
    fg.name = fascicleFileName;
end
return
