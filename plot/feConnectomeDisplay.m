function [fig_handle, h, fibers] = feConnectomeDisplay(varargin)
%
%	function [fig_handle, fibers, bundle_colors] = feConnectomeDisplay(
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
%     arg[2]= figure handle. Not optional. Should hte the handel to a
%             current figure.
%     arg[3]= (optional) A Nx3 matrix containing the colours in which the different bundles will be
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
%     return[1]= *fibers*, A 4xN cell structure containing the
%                 fiber data (to be used as arg[1] in subsequent
%                 calls of the function)
%
% 		return[2]= *bundle_colors*, A Nx3 matrix containing the colours of the
% 		            strands displayed (to be used as arg[3]
% 		            in subsequent calls of the function to keep colours
% 		            consistent). Colors are randomly generated if passed
% 		            empty.
%
%     A *fibers* is organized as follows:
%           fibers{i_fiber,2}  = index of the fiber in the collection.
%           fibers{i_fiber, 3} = the radius of the fiber proportion
%                                          of voxel.
%           fibers{i_fiber,4}  = bundle of fibers the fiber is part of.
%
%	USAGE:
%
%      feConnectomeDisplay(fg,figure);
%
%
% Franco (c) 2012 Stanford VISTA Team.
%
%	Based on the function 'dislay_strands' by Tom Close (2008) which was based on the 'tubeplot.m'
% by Anders Sandberg, asa@nada.kth.se, 2005).
tic
%% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

% Sort out inputs received and set default parameters.
[fibers,fig_handle, ~, ~, ~, ~,tube_corners, fiber_color] = checkInputs(varargin,nargin);
fprintf('\n[%s] Working on displaying %i fibers... ',mfilename, size(fibers,1));

minNodesNum = 3;

% Open the figure.
hold on;
box off;axis off;
%whitebg(fig_handle,[0 0 0]);
set(gcf,'Color',[0 0 0]);
set(gca,'color',[0 0 0])
num_fibers = size(fibers,1);  
%fiber_color = [.88 .5 .2];[.88 .9 .97];
X = cell(num_fibers,1);Y = X; Z = X;  segs = X;
fiberRad = .14135;
numNodes = zeros(num_fibers,1);

parfor i_fiber = 1:num_fibers
  % Get the current fiber and color
  fiber       = fibers{i_fiber,1};
  
  % Get the number of nodes, it is made of,
  % the number of segments the nodes define,
  % the length of each segment
  numNodes(i_fiber)    = size(fiber,1);
  segs{i_fiber}        = fiber(2: numNodes(i_fiber),:) - fiber(1: numNodes(i_fiber)-1,:);  
end

% Find the fibers that have more that a min number fo nodes
% We will plot only those
showme = find(numNodes >= minNodesNum);
t = cell(length(showme));n=t;b=t;
parfor i_fiber = 1:length(showme)
  % Calculate the frame for tube representing the fiber.
  [t{i_fiber},n{i_fiber},b{i_fiber}] = build3Dframe(fibers{showme(i_fiber)});

  % Build x,y,z coordinates for each node in the fibers and the
  % angle between them.
  [X{i_fiber},Y{i_fiber},Z{i_fiber}] = buildSurfaceFromFrame(fibers{showme(i_fiber)},fiberRad, ...
                                                             numNodes(showme(i_fiber)), ...
                                                             tube_corners, ...
                                                             segs{showme(i_fiber)}, ...
                                                             t{i_fiber},n{i_fiber},b{i_fiber});
end
clear t n b
h=nan(length(X),1);
for i_fiber = 1:length(X)
  % Make a surface for the fiber.
  h(i_fiber) = surf(X{i_fiber},Y{i_fiber},Z{i_fiber},'facecolor', fiber_color,'edgecolor', 'none');
end
% Format the figure.
formatFigure(fig_handle);

% Done
t = toc;
fprintf('done in  %2.3f seconds at %2.3f ms/fiber.\n',t,(t/num_fibers)*1000);

end % end MAIN function


%%%%%%%%%%%%%%%%
% formatFigure %
%%%%%%%%%%%%%%%%
function formatFigure(fig_handle)
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



%%%%%%%%%%%%%%%%%%%%%%%%%
% buildSurfaceFromFrame %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Z] = buildSurfaceFromFrame(fiber,fiber_radius,numNodes,tube_corners,segs,t,n,b)
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
  X(i_node,:) = fiber(i_node, 1) + (fiber_radius) * ( n_prime(1,1) * cos(theta) + b_prime(1,1) * sin(theta));
  Y(i_node,:) = fiber(i_node, 2) + (fiber_radius) * ( n_prime(1,2) * cos(theta) + b_prime(1,2) * sin(theta));
  Z(i_node,:) = fiber(i_node, 3) + (fiber_radius) * ( n_prime(1,3) * cos(theta) + b_prime(1,3) * sin(theta));
end

end % end function


%%%%%%%%%%%%%%%%
% build3Dframe %
%%%%%%%%%%%%%%%%
function [t,n,b] = build3Dframe(fiber)
%
% Build a tube frame that can be used for displaying.
%
% Franco
ref_vector = rand(1,3);
while (1)
  [t,n,b] = frame( fiber(:,1),fiber(:,2),fiber(:,3), ref_vector);
  if (all(~isnan(n))), break; end
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
  x=x'; y=y'; z=z';
  N=size(x,1);
end
p=[x y z];

t=zeros(N,3);b=zeros(N,3);n=zeros(N,3);
for i=2:(N-1)
  t(i,:)=(p(i+1,:)-p(i-1,:));
  tl=norm(t(i,:));
  if (tl>0), t(i,:)=t(i,:)/tl;
  else       t(i,:)=t(i-1,:);
  end
end

t(1,:)=p(2,:)-p(1,:);
t(1,:)=t(1,:)/norm(t(1,:));

t(N,:)=p(N,:)-p(N-1,:);
t(N,:)=t(N,:)/norm(t(N,:));

for i=2:(N-1)
  n(i,:)=cross(t(i,:),vec);
  nl=norm(n(i,:));
  if (nl>0), n(i,:)=n(i,:)/nl;
  else       n(i,:)=n(i-1,:);
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

end

%%%%%%%%%%%%%%%
% checkInputs %
%%%%%%%%%%%%%%%
function [fibers, fig_handle, collection_name, plot_indices, plot_3D,  plot_indices_type, ...
          tube_corners, bundle_colors] = checkInputs(var,nvarsin)
%
%
% Preprocess the inputs tot he main function and return the fibers and all the default
% parameters.
%

% Check the number of inputs of the main function.
if (nvarsin < 1 || nvarsin > 8)
  error('[%s] Incorrect number of arguments (supplied: %i)\n',mfilename,nvarsin);
end
radius = 0.2;
% Reorganize the fibers.
%
% Strip off first and last points from strand
numFibers = length(var{1}.fibers);
fibers    = cell(numFibers,1);
parfor ll = 1:numFibers
  fibers{ll,1} = var{1}.fibers{ll}';
  %fibers{ll,2} = ll - 1;
  %fibers{ll,4} = ll - 1;
  %fibers{ll,3} = radius;
end
collection_name = var{1}.name;

% Figure handle passed in by fePlot
if (nvarsin < 2), fig_handle = figure;
else   fig_handle = var{2};
end

% Colors for each bundle of fibers.
if (nvarsin < 3),  bundle_colors = [.88 .9 .97];
else bundle_colors = var{3};
end

% Indexes of the fibers that make up a bundle and are plotted in the same color.
if (nvarsin >= 4), plot_indices = var{4};
else plot_indices = [];
end

% Name for the group of fibers plotted together, default is bundle,
% but it can also be the name of a fiber group.
if (nvarsin >= 5), plot_indices_type = var{5};
else plot_indices_type = 'bundle';
end

% Choose whether to plot simple lines or 3D lines.
if (nvarsin >= 6 && strcmp('line', var{6})), plot_3D = 0;
else plot_3D = 1;
end

% Choose the resolution at which the fiber section s samples
% Low when many fibers are passed in, high otherwise.
if (nvarsin >= 7), tube_divs = var{7};
else
  if (size(fibers,1) <= 100),     tube_divs = 160;
  elseif (size(fibers,1) <= 500), tube_divs = 20;
  else tube_divs = 10;
  end
end

% the number of corners in eahc tube.
tube_corners = tube_divs + 1;
bundle_indices = [];

% Make colors for each bundle.= of strands.
% if (size(bundle_colors,1) == 0)
%   %   parfor i_fiber = 1:size(fibers,1)
%   %     bundle_indices = [bundle_indices; fibers{i_fiber, 4}];
%   %   end
%   bundle_indices = unique(bundle_indices);
%   bundle_indices = sort(bundle_indices);
%   bundle_colors  = rand([(max(bundle_indices)+1), 3]);
% end

end % end function


