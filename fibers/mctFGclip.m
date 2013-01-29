function fg = mctFGclip(fg,coords)
%
% function fg = mctFGclip(fg,coords)
%
% Clips fibers in a fiber group to the coordinates in coords.
% 
% Ideally the coordinates are those of an roi, the operation will limit the
% fiber paths to within the ROI.
%
% fg and coords are assumed to be in the same coordinate space (acpc or
% image).
%
% See also:
%
% Example:
%
% Franco & Ariel
%
% (c) Stanford VISTA Team 2012

fprintf('\n[%s] Clipping fibers to be constrained within the volume\n        defined by the passed coordinates.\n',mfilename)

% Get the nodes to voxel matching. This is a cell array size(fg.fibers),
% each entry in the cell array tells you which node is in which voxel.
% 
nodes2voxels = fefgGet(fg,'nodes2voxels',coords);

% Find the indexes of the nodes that do pass through the roi
% coordinates.
nodesIndex = cellfun(@find,nodes2voxels, 'UniformOutput', false);

% Check whether the indexes of the nodes are consecutives. If they are not
% it means that they went out of the ROI then back in the ROI.
gradient = cellfun(@diff,nodesIndex, 'UniformOutput', false);

% If a fiber only has one node in the ROI gradient will be emtpy for that
% fiber. We set the fiber to 1.
gradient = cellfun(@oneifempty,gradient, 'UniformOutput', false);

% Now we deal with cases in which the fibers go out of the roi and then go
% back in into the ROI in these cases we split the fibers into separated
% fibers.
c = 1; % THis is the counter that keeps track of the new fibers. 
       % The fiber number will increase after the following operations
for iif = 1:length(fg.fibers)
  % Find the nodes that are continous within the ROI
  isone = gradient{iif} == 1;
  
  % If all the nodes in this fiber were in the ROI just copy the fiber
  if all(isone)
    fibers{c} = fg.fibers{iif}(:,nodesIndex{iif});
    if isfield(fg,'seeds') && ~isempty(fg.seeds)
      if size(fg.seeds,1) == 3
        seeds(1:3,c)  = fg.seeds(:,iif);
      elseif size(fg.seeds,2) == 3
        seeds(c,1:3)  = fg.seeds(iif,:);
      end
    end
    if isfield(fg,'Q')
      Q{c}  = fg.Q{iif};
    end
    
    c = c + 1;
    
  % If the fiber exits the ROI and then re-enters it, ctu the fibers into
  % segements and save the segments as separate fibers.
  else
    indx = find(~isone);
    for iin = 1:length(indx)
      if iin == 1
        fibers{c} = fg.fibers{iif}(:, nodesIndex{iif}(1:indx(iin)));
        if isfield(fg,'seeds') && ~isempty(fg.seeds)
          seeds(1:3,c)  = fg.seeds(:,iif);
        end
        if isfield(fg,'Q')
          Q{c}  = fg.Q{iif};
        end
      else
        fibers{c} = fg.fibers{iif}(:, nodesIndex{iif}(indx(iin - 1) + 1:indx(iin)));
        if isfield(fg,'seeds') && ~isempty(fg.seeds)
          seeds(1:3,c)  = fg.seeds(:,iif);
        end
        if isfield(fg,'Q')
          Q{c}  = fg.Q{iif};
        end
      end
      c = c + 1;
    end
  end
end

% Save to output fiber group
fg.name = sprintf('%s clipped to ROI coords',fg.name);
if isfield(fg,'Q')
  fg.Q = Q;
end
if isfield(fg,'seeds') && ~isempty(fg.seeds)
  fg.seeds = seeds;
end

fg.fibers = fibers;

% The next line is just fo rdebugging if the code is correct we should
% never have a problme withthe following line.
% Delete me soonish.
deleteme = fefgGet(fg,'nodes2voxels',coords);



%-------------------------%
function x = oneifempty(x)
%
% This function is used to substitute empty matrices in a cell with 1.
%
if isempty(x), x = 1;end