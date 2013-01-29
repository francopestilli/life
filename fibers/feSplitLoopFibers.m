function fg = feSplitLoopFibers(fg)
%
% function fg = feSplitLoopFibers(fg)
%
% Splits fibers in a fiber group that are discontinuous.
% 
% These fibers are probably entering and existing a volume ROI. This type
% of fibers are often the result of a a clipping process (feClipFibersTovolume.m)
%
% See also: feClipFibersToVolume.m
%
% Example:
%
% Franco (c) Stanford VISTA Team 2012

fprintf('\n[%s] Splitting fibers that are discontinuous.\n',mfilename)

% Check whether the nodes are consecutives. If they are not
% it means that they went out of the ROI then back in the ROI.
% Compute the distance between nodes, this is generally a small number (0.5).
% But for fibers that exitedna dnre-entered the ROI it will be large.
dx = cellfun(@getNodeDistance,fg.fibers, 'UniformOutput', false);

% Compute the median distance. This is rounds up and then used as a
% thresholds for detecting the fibers and the nodes that exited and
% reentered.
md = floor(100*median(cell2mat(dx')))/100;

% Now we arbitrarily decide that we keep spli only the fibers that have
% nodes with double this median distnce, this is because if we have fibers
% that exit and reenter quickly they will show up OK on display.
% This small hack speeds up the code and the visualization si faster
maxDistance = 2*md - 0.02;

% If a fiber only has one node in the ROI dx will be emtpy for that
% fiber. We set the fiber to 1.
dx = cellfun(@oneifempty,dx, 'UniformOutput', false);

% Precompute the final number of fibers, the ones we had plus the new ones
% that we will create by splitting the fibers that are discontinuos.
c = 1;
parfor iif = 1:length(fg.fibers)
  % Find the nodes that are continous within the ROI
  isone = dx{iif} <= maxDistance;
  if all(isone), c = c + 1;
  else           c = c + length(find(~isone));
  end
end
% Use the total number fo fibers to initialize the fibers and the other
% fields
fibers = cell(c,1); Q = fibers;
if     (size(fg.seeds,1) == 3), seeds = nan(3,c);
elseif (size(fg.seeds,2) == 3), seeds = nan(c,3);
end

% Now we deal with cases in which the fibers go out of the roi and then go
% back in into the ROI in these cases we split the fibers into separated
% fibers.
c = 1; % This is the counter that keeps track of the new fibers. 
       % The fiber number will increase after the following operations
for iif = 1:length(fg.fibers)
  % Find the nodes that are continous within the ROI
  isone = dx{iif} <= maxDistance;
  
  % If all the nodes in this fiber were in the ROI just copy the fiber
  if all(isone)
    fibers{c} = fg.fibers{iif};
%     if isfield(fg,'seeds') && ~isempty(fg.seeds)
%       if size(fg.seeds,1) == 3
%         seeds(1:3,c)  = fg.seeds(:,iif);
%       elseif size(fg.seeds,2) == 3
%         seeds(c,1:3)  = fg.seeds(iif,:);
%       end
%     end
%     if isfield(fg,'Q')
%       Q{c}  = fg.Q{iif};
%     end
    
    c = c + 1;
    
  % If the fiber exits the ROI and then re-enters it, ctu the fibers into
  % segements and save the segments as separate fibers.
  else
    indx = find(~isone);
    for iin = 1:length(indx)
      if iin == 1
        fibers{c} = fg.fibers{iif}(:, (1:indx(iin)));
%         if isfield(fg,'seeds') && ~isempty(fg.seeds)
%           seeds(1:3,c)  = fg.seeds(:,iif);
%         end
%         if isfield(fg,'Q')
%           Q{c}  = fg.Q{iif};
%         end
      else
        fibers{c} = fg.fibers{iif}(:, (indx(iin - 1) + 1:indx(iin)));
%         if isfield(fg,'seeds') && ~isempty(fg.seeds)
%           seeds(1:3,c)  = fg.seeds(:,iif);
%         end
%         if isfield(fg,'Q')
%           Q{c}  = fg.Q{iif};
%         end
      end
      c = c + 1;
    end
  end
end

% Save to output fiber group
fg.name = sprintf('%s clipped to ROI coords',fg.name);
% if isfield(fg,'Q')
%   fg.Q = Q;
% end
% if isfield(fg,'seeds') && ~isempty(fg.seeds)
%   fg.seeds = seeds;
% end
fg.fibers = fibers;

%-------------------------%
function stepSize = getNodeDistance(fiber)
% Computes the distance between nodes in a ifber.
stepSize = sqrt(sum(diff(fiber,1,2).^2));


%-------------------------%
function x = oneifempty(x)
%
% This function is used to substitute empty matrices in a cell with 1.
%
if isempty(x), x = 1;end