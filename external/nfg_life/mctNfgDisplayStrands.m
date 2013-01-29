function [fig_handle, fiber_collection, bundle_colors] = mctNfgDisplayStrands(varargin)
%
%	function [fig_handle, fiber_collection, bundle_colors] = mctNfgDisplayStrands(
%                                                           fiber group structure, ...
%                                                           strands dirname, or strand structure,
%                                                           [bundle_colors],
%                                                           [bundle or strands to plot],
%                                                           ['bundle' or 'fiber' (which ones to plot)],
%                                                           ['3D' or 'line' (for 3D rendered or line plots)],
%                                                           [tube_divs, resolution at which the tube is rendered].
%                                                           [radius = size of the fibers, either one numbe or a vector 
%                                                                     of the size of fibers]);
%	
% INPUTS:
%     arg[1]= Either:
%             1. the directory name of the stored collection
%             2. or a cell structure containing the loaded strands from a previous call (return variable 1).
%             3. or a fiber group.
%
%     arg[2]= (optional) A Nx3 matrix containing the colours in which the different bundles will be
% 		        displayed.
% 	  arg[4]= (optional) The indices of the strands/bundles to be displayed.
% 	  arg[5]= (optional) Whether the supplied indices correspond to 'fiber' or 'bundle' (default
% 		        'bundle')
% 	  arg[6]= (optional) The number of divisions around the diameter of the fiber. The default is
%             6, which corresponds to a hexagonal cross-section.
% 	  arg[7]= (optional) Radius of each fiber in a fiber group (this will only work if a fiber group is passed in)
%
% OUTPUTS:
%
%     return[1]= *fiber_collection*, A 4xN cell structure containing the
%                 fiber data (to be used as arg[1] in subsequent
%                 calls of the function)
%
% 		return[2]= *bundle_colors*, A Nx3 matrix containing the colours of the
% 		            strands displayed (to be used as arg[3]
% 		            in subsequent calls of the function to keep colours
% 		            consistent). Colors are randomly generated if passed
% 		            empty.
%
%     A *fiber_collection* is organized as follows:
%           fiber_collection{i_fiber,2}  = index of the fiber in the collection.
%           fiber_collection{i_fiber, 3} = the radius of the fiber proportion
%                                          of voxel.
%           fiber_collection{i_fiber,4}  = bundle of fibers the fiber is part of.
%
%	USAGE:
%
%      mctNfgDisplayStrands(fg);
%
%
% Franco Pestilli
%
%	Based on the function 'dislay_strands' by Tom Close (2008) which was based on the 'tubeplot.m' 
% by Anders Sandberg, asa@nada.kth.se, 2005).

global bundle_colors;  % bundle colours is made global in case you forget to save the current colours and want to retrieve them.

minNumOfNodes = 2; % This sets the minimum number of nodes for afiber to be displayed.

% Sort out inputs received and set default parameters.
[fiber_collection,   ...
  collection_name,   ...
  isotropic_regions, ...
  plot_indices,      ...
  plot_3D,           ...
  plot_indices_type, ...
  tube_corners,      ...
  bundle_colors] = checkInputs(varargin,nargin);

fprintf('\n[%s] Working on displaying %i strands.\n',mfilename, size(fiber_collection,1));

% Open the figure.
fig_handle = mrvNewGraphWin(collection_name); hold on;
box off;axis off;
whitebg(fig_handle,[0 0 0]); 
set(fig_handle,'Color',[0 0 0])

% Plot isotropic regions if any.
if ~isempty(isotropic_regions)
  makeIsotropicSurface(isotropic_regions);
end

% Build a tube as 3D surface (or just plot aline) for each fiber.
fiber_counter = 0;
for i_fiber = 1:size(fiber_collection,1)
  if (size(fiber_collection{i_fiber,1},1) < minNumOfNodes)
   % fprintf('\n[%s] Strand %i is not shown because it only has %i control points.\n', ...
   %   mfilename, i_fiber,size(fiber_collection{i_fiber,1},1));
    
  elseif ( (strcmp(plot_indices_type,'bundle') && ...
      any(fiber_collection{i_fiber,4} == plot_indices)) || ...
      (strcmp(plot_indices_type,'strand') && ...
      any(fiber_collection{i_fiber,2} == plot_indices)) || ...
      size(plot_indices,1) == 0)
    
    % Get the current fiber and color
    fiber = fiber_collection{i_fiber,1};
    fiber_color = bundle_colors(fiber_collection{i_fiber, 4}+1,:);
    
    % Plot it.
    if (plot_3D)
      % get the fiber index, it radius and the index of the bundle it is part of
      save_fiber_index  = fiber_collection{i_fiber,2};
      fiber_radius      = fiber_collection{i_fiber, 3};
      save_bundle_index = fiber_collection{i_fiber,4};
      
      % Get the number of nodes, it is made of,
      % the number of segments the nodes define,
      % the length of each segment
      numNodes    = size(fiber,1);
      segs        = fiber(2:numNodes,:) - fiber(1:numNodes-1,:);
      length_segs = sqrt(dot(segs,segs,2));
      
      % Set the legnth of every segment that is zero-length to a very
      % small but plottable length.
      zero_length_segs = find(length_segs == 0.0);
      if (size(zero_length_segs,1) >= 1)
        length_segs(zero_length_segs) = 0.0001;
        segs(zero_length_segs) = 0.0001;
        fiber(zero_length_segs) = fiber(zero_length_segs)  + 0.0001;
      end
      
      % Calculate the frame for tube representing the fiber.
      [t,n,b] = buildFrameFromStrands(fiber);
      
      % Build x,y,z coordinates for each node in the fibers and the
      % angle between them.
      if length(fiber_radius) == size(fiber_collection,1)
        fiberRad = fiber_radius(i_fiber);
      else
        fiberRad = fiber_radius;
      end
      [X,Y,Z] = buildSurfaceCoordsFromStrandsFrame(fiber,fiberRad,numNodes,tube_corners,segs,t,n,b);
      
      % Make a surface for the fiber.
      surf(X,Y,Z,'facecolor', fiber_color,'edgecolor', 'none');
      
    else % If we are not plotting 3D simply plot lines.
      plot3(fiber(:,1), fiber(:,2), fiber(:,3), 'Color', fiber_color);
    end
    
    % Update the counter for the current fiber.
    fiber_counter = fiber_counter + 1;
    % Display the status of full process. one '.' each 25 strands.
    if (mod(fiber_counter,25) == 0), fprintf('.'); end
  end
end

% Add the name to the fiber group.
fiber_collection{end+1,1} = collection_name;

% Format the figure.
formatFigure(fig_handle,fiber_counter,size(fiber_collection,1));

% Done
fprintf('\nDone displaying %i strands. \n\n',fiber_counter);

end % end MAIN function


%%%%%%%%%%%%%%%%
% formatFigure %
%%%%%%%%%%%%%%%%
function formatFigure(fig_handle,fiberNum,totFiberNum)
%
% Format the figure by adding menus,
% lightning, lables, etc.
%

% Format figure.
hold off;
set(fig_handle,'Units','normalized', 'Position', [0 0.2 0.2 0.4], 'DoubleBuffer', 'on','Color',[0 0 0]);
cameratoolbar('Show');
cameratoolbar('SetMode','orbit');
light
lighting gouraud;

%title(sprintf('Displayed the %i longer fibers out of %i total fibers.',fiberNum,totFiberNum))

% Format axis
grid on;
label_handle = get(gca,'xlabel');
set(label_handle, 'string', 'X-axis', 'color', [0 0 0]);
label_handle = get(gca,'ylabel');
set(label_handle, 'string', 'Y-axis', 'color', [0 0 0]);
label_handle = get(gca,'zlabel');
set(label_handle, 'string', 'Z-axis', 'color', [0 0 0]);
daspect ([ 1 1 1 ]);

end % end function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% buildSurfaceCoordsFromStrandsFrame %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Z] = buildSurfaceCoordsFromStrandsFrame(fiber,fiber_radius,numNodes,tube_corners,segs,t,n,b)
%
% Build the coordinates of a 3D surface that can be used with surf.m,
% from the frames of each fiber.
%
% Franco
X = zeros(numNodes, tube_corners);
Y = zeros(numNodes, tube_corners);
Z = zeros(numNodes, tube_corners);

theta = 0 : (2*pi/(tube_corners-1)) : (2*pi);
for i_node = 1:numNodes
  if (i_node == 1)
    w       = fiber(1, :) + n(1,:);
    n_prime = n(1,:);
    b_prime = b(1,:);
  else
    mu      = dot(t(i_node,:), fiber(i_node,:) - w, 2) / dot(t(i_node,:), segs(i_node-1,:),2);
    w_proj  = w + mu * segs(i_node-1, :);
    n_prime = w_proj - fiber(i_node,:);
    n_prime = n_prime ./ norm(n_prime);
    b_prime = cross( t(i_node,:), n_prime);
    b_prime = b_prime ./ norm(b_prime);
    w       = fiber(i_node,:) + n_prime;
    
  end
  X(i_node,:) = fiber(i_node, 1) + (fiber_radius + randn(1,1)*.01) * ( n_prime(1,1) * cos(theta) + b_prime(1,1) * sin(theta));
  Y(i_node,:) = fiber(i_node, 2) + (fiber_radius + randn(1,1)*.01) * ( n_prime(1,2) * cos(theta) + b_prime(1,2) * sin(theta));
  Z(i_node,:) = fiber(i_node, 3) + (fiber_radius + randn(1,1)*.01) * ( n_prime(1,3) * cos(theta) + b_prime(1,3) * sin(theta));
end

end % end function


%%%%%%%%%%%%%%%%%%%%%%%%%
% buildFrameFromStrands %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,n,b] = buildFrameFromStrands(fiber)
%
% Build a tube frame that can be used for displaying.
%
% Franco
count = 0;
while (1)
  ref_vector = rand(1,3);
  [t,n,b] = frame( fiber(:,1),fiber(:,2),fiber(:,3), ref_vector);
  % Exit if any fiber as a NAN in it.
  if (all(~isnan(n))), break; end
  if (count > 100), error('likely nan values in fiber %s!',i_fiber-1);end
  count = count + 1;
end

end % end function


%%%%%%%%%%%%%%%
% checkInputs %
%%%%%%%%%%%%%%%
function [fiber_collection, ...
  collection_name,   ...
  isotropic_regions, ...
  plot_indices,      ...
  plot_3D,           ...
  plot_indices_type, ...
  tube_corners,         ...
  bundle_colors] = checkInputs(var,nvarsin)
%
%
% Preprocess the inputs tot he main function and return the fibers and all the default
% parameters.
%

% Check the number of inputs of the main function.
if (nvarsin < 1 || nvarsin > 7)
  error('[%s] Incorrect number of arguments (supplied: %i)\n',mfilename,nvarsin);
end

isotropic_regions = [];

% Check what kind of input we received for the fibers.
% Fibers are LiFe fibers.
% Strands are NFG fibers.
if (ischar(var{1})) % a directory to load from a strand collection
  %[fiber_collection, isotropic_regions] = load_strand_collection(var{1});
  [fiber_collection, isotropic_regions] = mctNfgLoadStrands(var{1});
  collection_name = var{1};
  
elseif (iscell(var{1})) % a strand collection
  
  fiber_collection = var{1};
  if size(fiber_collection,1) == 0
    error('No strands in loaded collection');
  end
  collection_name = fiber_collection{end,1};
  fiber_collection(end,:) = [];
  
elseif isstruct(var{1}) % fiber group
  
  if (nvarsin >= 7), radius = (var{7} + 0.06) * 2; % here I add a little bit so 
                                                   % that 0-thickness fibers are actually 
                                                   % visible and I scale up
                                                   % the passed thickness
                                                   % so that fiber weights
                                                   % used as thickness are
                                                   % visible
  else radius = 0.5;
  end
  
  % transform the fiber group into a fiber collection
  fiber_collection = mctNfgFibers2Strands(var{1},radius);
  collection_name = fiber_collection{end,1};
  fiber_collection(end,:) = [];
  
else
  error('Incorrect data type for argument 1\n');
end

% Colors for each bundle of fibers.
if (nvarsin >= 2),  bundle_colors = repmat(var{2},size(fiber_collection,1),1);
else bundle_colors = [];
end

% Indexes of the fibers that make up a bundle and are plotted in the same color.
if (nvarsin >= 3), plot_indices = var{3};
else plot_indices = [];
end

% Name for the group of fibers plotted together, default is bundle,
% but it can also be the name of a fiber group.
if (nvarsin >= 4), plot_indices_type = var{4};
else plot_indices_type = 'bundle';
end

% Choose whether to plot simple lines or 3D lines.
if (nvarsin >= 5 && strcmp('line', var{5})), plot_3D = 0;
else plot_3D = 1;
end

% Choose the resolution at which the fiber section s samples
% Low when many fibers are passed in, high otherwise.
if (nvarsin >= 6)
  tube_divs = var{6};
else
  if (size(fiber_collection,1) <= 100)
    tube_divs = 160;
  elseif (size(fiber_collection,1) <= 300)
    tube_divs = 60;
  else tube_divs = 40;
  end
end

% the number of corners in eahc tube.
tube_corners = tube_divs + 1;

display_key = (size(bundle_colors,1) == 0);
bundle_indices = [];

% Make colors for each bundle.= of strands.
if (size(bundle_colors,1) == 0)
  for i_fiber = 1:size(fiber_collection,1)
    bundle_indices = [bundle_indices; fiber_collection{i_fiber, 4}]; %#ok<AGROW>
  end
  bundle_indices = unique(bundle_indices);
  bundle_indices = sort(bundle_indices);
  bundle_colors  = rand([(max(bundle_indices)+1), 3]);
end

% Display a figure with the fiber indexes.
if (display_key)
  %display_colour_key(bundle_colors, bundle_indices);
end

end % end function


%%%%%%%%%
% frame %
%%%%%%%%%
function [t,n,b] = frame(x,y,z,vec)
%
% FRAME Calculate a Frenet-like frame for a polygonal space curve
% [t,n,b]=frame(x,y,z,v) returns the tangent unit vector, a normal
% and a binormal of the space curve x,y,z. The curve may be a row or
% column vector, the frame vectors are each row vectors.
%
% This function calculates the normal by taking the cross product
% of the tangent with the vector v; if v is chosen so that it is
% always far from t the frame will not twist unduly.
%
% If two points coincide, the previous tangent and normal will be used.
%
% Written by Anders Sandberg, asa@nada.kth.se, 2005



N=size(x,1);
if (N==1)
  x=x'; %'
  y=y'; %'
  z=z'; %'
  N=size(x,1);
end

t=zeros(N,3);
b=zeros(N,3);
n=zeros(N,3);

p=[x y z];

for i=2:(N-1)
  t(i,:)=(p(i+1,:)-p(i-1,:));
  tl=norm(t(i,:));
  if (tl>0)
    t(i,:)=t(i,:)/tl;
  else
    t(i,:)=t(i-1,:);
  end
end

t(1,:)=p(2,:)-p(1,:);
t(1,:)=t(1,:)/norm(t(1,:));

t(N,:)=p(N,:)-p(N-1,:);
t(N,:)=t(N,:)/norm(t(N,:));

for i=2:(N-1)
  n(i,:)=cross(t(i,:),vec);
  nl=norm(n(i,:));
  if (nl>0)
    n(i,:)=n(i,:)/nl;
  else
    n(i,:)=n(i-1,:);
  end
end

n(1,:)=cross(t(1,:),vec);
n(1,:)=n(1,:)/norm(n(1,:));

n(N,:)=cross(t(N,:),vec);
n(N,:)=n(N,:)/norm(n(N,:));

for i=1:N
  b(i,:)=cross(t(i,:),n(i,:));
  b(i,:)=b(i,:)/norm(b(i,:));
end

end % end function


%%%%%%%%%%%%%%%%%%%%%%%%
% makeIsotropicSurface %
%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeIsotropicSurface(isotropic_regions)
% Build a surface for each isotropic Region.
% Isoptropic regions are returned as part fo the NFG simuation.
% The simulate CSF.
for isotropic_region_i = 1 : size(isotropic_regions,1)
  centre = isotropic_regions(isotropic_region_i, 1:3)';
  radius = isotropic_regions(isotropic_region_i, 4);
  num_faces = round(radius * 100);
  [X, Y, Z] = sphere(num_faces);
  X = X * radius + centre(1);
  Y = Y * radius + centre(2);
  Z = Z * radius + centre(3);
  
  % Surfi it out...
  h = surf(X, Y, Z);
  set(h,'facecolor', [0.125 0.125 0.125]);
  set(h, 'FaceAlpha', 0.125);
end

end % end function


