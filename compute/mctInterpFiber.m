function intFiber = mctInterpFiber(res,fiber)
% fiber should be numNodesX3

% computer the fiber length.
len = sum(sqrt(sum(diff(fiber,[],1).^2,2)));

% compute the number of nodes in this fiber at this resolution given the
% fiber length
numNodes = ceil(len*res);
%fprintf('\nFiber Length: %2.2f and num Nodes: %i, ratio: %2.4f\n',len, numNodes,len/numNodes)

numNodes = linspace(0,1,numNodes)';

% how many nodes will be interpolated?
intNumNodes = numel(numNodes);

% the number of nodes in the fibers passed in.
% assumes the fiber is 3xn, where the first dimension contains the x,y,z
% coordinates in the 3D volume.
curNumNodes = size(fiber,1);

% preallocate the result, fiber
intFiber = NaN(intNumNodes,3);

% Compute the chordal (linear) arclength
% of each segment. This will be needed for
% any of the methods.
chordlen = sqrt(sum(diff(fiber,[],1).^2,2));

% Normalize the arclengths to a unit total
chordlen = chordlen/sum(chordlen);

% cumulative arclength
cumarc = [0;cumsum(chordlen)];

% compute parametric splines
spl = cell(1,3);
spld = spl;
diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
for i = 1:3
  spl{i} = pchip(cumarc,fiber(:,i));
  
  % and now differentiate them
  xp = spl{i};
  xp.coefs = xp.coefs*diffarray;
  xp.order = 3;
  spld{i} = xp;
end

% Generate the total arclength along the curve
% by integrating each segment and summing the
% results. The integration scheme does its job
% using an ode solver.

% polyarray here contains the derivative polynomials
% for each spline in a given segment
polyarray = zeros(3,3);
seglen = zeros(curNumNodes-1,1);

% options for ode45
opts = odeset('reltol',1.e-9);
for i = 1:spl{1}.pieces
  % extract polynomials for the derivatives
  for j = 1:3
    polyarray(j,:) = spld{j}.coefs(i,:);
  end
  
  % integrate the arclength for the i'th segment
  % using ode45 for the integral. I could have
  % done this part with quad too, but then it
  % would not have been perfectly (numerically)
  % consistent with the next operation in this tool.
  [tout,yout] = ode45(@(numNodes,y) segkernel(numNodes,y),[0,chordlen(i)],0,opts); %#ok
  seglen(i) = yout(end);
end

% and normalize the segments to have unit total length
totalsplinelength = sum(seglen);
cumseglen = [0;cumsum(seglen)];

% which interval did each point fall into, in
% terms of numNodes, but relative to the cumulative
% arc lengths along the parametric spline?
[~,tbins] = histc(numNodes*totalsplinelength,cumseglen);

% catch any problems at the ends
tbins((tbins <= 0) | (numNodes <= 0)) = 1;
tbins((tbins >= curNumNodes) | (numNodes >= 1)) = curNumNodes - 1;

% Do the fractional integration within each segment
% for the interpolated points. numNodes is the parameter
% used to define the splines. It is defined in terms
% of a linear chordal arclength. This works nicely when
% a linear piecewise interpolant was used. However,
% what is asked for is an arclength interpolation
% in terms of arclength of the spline itself. Call s
% the arclength traveled along the spline.
s = totalsplinelength*numNodes;

% the ode45 options will now include an events property
% so we can catch zero crossings.
opts = odeset('reltol',1.e-9,'events',@ode_events);
for i = 1:intNumNodes
  % si is the piece of arc length that we will look
  % for in this spline segment.
  si = s(i) - cumseglen(tbins(i));
  
  % extract polynomials for the derivatives
  % in the interval the point lies in
  for j = 1:3
    polyarray(j,:) = spld{j}.coefs(tbins(i),:);
  end
  
  % we need to integrate in numNodes, until the integral
  % crosses the specified value of si. Because we
  % have defined totalsplinelength, the lengths will
  % be normalized at this point to a unit length.
  %
  % Start the ode solver at -si, so we will just
  % look for an event where y crosses zero.
  [tout,yout,te,ye] = ode45(@(numNodes,y) segkernel(numNodes,y),[0,chordlen(tbins(i))],-si,opts); %#ok
  
  % we only need that point where a zero crossing occurred
  % if no crossing was found, then we can look at each end.
  if ~isempty(te)
      ti = te(1) + cumarc(tbins(i));
  else
    % a crossing must have happened at the very
    % beginning or the end, and the ode solver
    % missed it, not trapping that event.
    if abs(yout(1)) < abs(yout(end))
      % the event must have been at the start.
      ti = tout(1) + cumarc(tbins(i));
    else
      % the event must have been at the end.
      ti = tout(end) + cumarc(tbins(i));
    end
  end
  
  % Interpolate the parametric splines at ti to get
  % our interpolated value.
  for j = 1:3
    intFiber(i,j) = ppval(spl{j},ti);
  end
  
end


% ===============================================
%  nested function for the integration kernel
% ===============================================
  function val = segkernel(numNodes,y) %#ok
    % sqrt((dx/dt)^2 + (dy/dt)^2 + ...)
    val = zeros(size(numNodes));
    for k = 1:3
      val = val + polyval(polyarray(k,:),numNodes).^2;
    end
    val = sqrt(val);
    
  end % function segkernel

% ===============================================
%  nested function for ode45 integration events
% ===============================================
  function [value,isterminal,direction] = ode_events(numNodes,y) %#ok
    % ode event trap, looking for zero crossings of y.
    value = y;
    isterminal = ones(size(y));
    direction = ones(size(y));
  end % function ode_events

end % Main function - mctInterpFiber


