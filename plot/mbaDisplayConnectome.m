function [figureHandle, lightHandle] = mbaDisplayConnectome(fibers,figureHandle,fiberColor,colorType,cMap, fiberRadius,numSurfaceFaces)
%
%	function [figureHandle, lightHandle] = mbaDisplayConnectome(fibers,figureHandle,color,fiberRadius,numSurfaceFaces);
%
% INPUTS:
%     fibers         - A fiber group.
%     fiureHandle - figure handle. Optional.
%     fiberColor  - A Nx3 matrix containing the colours in which the different bundles will be
% 		        displayed.
% 	  colorType   - 'map' or 'single'. Optional.
% 	  cMap        - A matlab colormap (e.g., map=hsv)
% 	  fiberRadius - The thicknes of each fiber.
%     numSurfaceFaces - The number of divisions around the diameter of the fiber.
%                       Default = 10.
%
% OUTPUTS:
%     figureHandel - Handle to the figure with the Connectome.
% 	   lightHandle - HAndle to the current ligth used on the connectome.
%
%	USAGE:
%      mbaDisplayConnectome(fg,figure);
%
% Written by Franco Pestilli (c) Stanford University, 2013.

tic
%% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

% Check if a fgure handle was passed in, otherwise open a figure
if notDefined('figureHandle'), figureHandle = figure;end

% Choose whether to plot simple lines or 3D lines.
if notDefined('plot2d') plot2d = 0;end

% Choose the resolution at which the fiber section s samples
% Low when many fibers are passed in, high otherwise.
if notDefined('numSurfaceFaces'), numSurfaceFaces = 50; end
if notDefined('fiberColor'),  
    fiberColor = [.84,.83,.99];
    colorType = 'single';
end
if notDefined('fiberRadius'),fiberRadius = .5;end
if notDefined('minNodesNum'), minNodesNum = 3; end

% This is the number of edges eac surface has
surfaceCorners = numSurfaceFaces + 1;

% Reorganize the fibers.
if (size(fibers{1},1) == 3) 
   fibers =  cellfun(@transpose,fibers,'UniformOutput',0);
end
numFibers = length(fibers);
X = cell(numFibers,1);
Y = X; Z = X;  segs = X;
numNodes = zeros(numFibers,1);

% We remove the first and last nodes from each fiber and colelct the number
% of nodes per fiber:
parfor i_fiber = 1:numFibers
 
  % Get the number of nodes, it is made of,
  % the number of segments the nodes define,
  % the length of each segment
  numNodes(i_fiber) = size(fibers{i_fiber},1);
  segs{i_fiber}     = fibers{i_fiber}(2: numNodes(i_fiber),:) - fibers{i_fiber}(1: numNodes(i_fiber)-1,:);  
end

% Find the fibers that have more that a min number fo nodes
% We will plot only those
showme  = find(numNodes >= minNodesNum);

%t = cell(length(showme));n=t;b=t;
parfor i_fiber = 1:length(showme)
    % Make a variable fiber radius:
    fr = fiberRadius + (fiberRadius/10) .* randn(size(fibers{showme(i_fiber)}));
    
    % Calculate the frame for tube representing the fiber.
    [t,n,b] = build3Dframe(fibers{showme(i_fiber)});
    
    % Build x,y,z coordinates for each node in the fibers and the
    % angle between them.
    [X{i_fiber},Y{i_fiber},Z{i_fiber}] =                   ...
        buildSurfaceFromFrame(fibers{showme(i_fiber)},     ...
        fr, ...
        numNodes(showme(i_fiber)),   ...
        surfaceCorners,              ...
        segs{showme(i_fiber)},       ...
        t,n,b);
end
clear t n b

% Now if we received only one color we generate a cell array with the same
% size of the fibers by replicating the color for each fiber/node:
h=nan(length(X),1);  hold on
switch colorType
    case {'single','uniform'}
        for i_fiber = 1:length(X)
            % Make a surface for the fiber.
            h(i_fiber) = surf(X{i_fiber},Y{i_fiber},Z{i_fiber}, 'facecolor',fiberColor,'edgecolor', 'none');
        end
 
    case {'map'}
        for i_fiber = 1:length(showme)
            % Make a surface for the fiber.
            h(i_fiber) = surf(X{i_fiber},Y{i_fiber},Z{i_fiber}, repmat(fiberColor{showme(i_fiber)},1,size(Z{i_fiber},2)), ...
                'edgecolor', 'none',  'AmbientStrength',0.9);
        end
        colormap(cMap)
    otherwise
        keyboard
end
    
% Format the figure.
lightHandle = formatFigure(gcf);

% Done
t = toc;
fprintf('[%s] done in  %2.3f seconds at %2.3f ms/fiber.\n',mfilename, t,(t/numFibers)*1000);

end % end MAIN function


%%%%%%%%%%%%%%%%
% formatFigure %
%%%%%%%%%%%%%%%%
function lh = formatFigure(fig_handle)
%
% Format the figure by adding menus,
% lightning, lables, etc.
%

% Format figure.
box off;axis off;
set(gcf,'Color',[0 0 0]);
set(gca,'color',[0 0 0])
hold off;
set(fig_handle, ...
    'Units','normalized', ...
    'Position', [0 0.2 0.2 0.4], ...
    'DoubleBuffer', 'on', ...
    'Color',[0 0 0]);

cameratoolbar('Show');
cameratoolbar('SetMode','orbit');
lh = light;
lighting gouraud;

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
function [X,Y,Z] = buildSurfaceFromFrame(fiber,fiber_radius,numNodes,surfaceCorners,segs,t,n,b)
%
% Build the coordinates of a 3D surface that can be used with surf.m,
% from the frames of each fiber.
%
% Franco
X = zeros(numNodes, surfaceCorners);
Y = zeros(numNodes, surfaceCorners);
Z = zeros(numNodes, surfaceCorners);

theta = 0 : (2*pi/(surfaceCorners-1)) : (2*pi);
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
  X(i_node,:) = fiber(i_node, 1) + (fiber_radius(i_node)) * ( n_prime(1,1) * cos(theta) + b_prime(1,1) * sin(theta));
  Y(i_node,:) = fiber(i_node, 2) + (fiber_radius(i_node)) * ( n_prime(1,2) * cos(theta) + b_prime(1,2) * sin(theta));
  Z(i_node,:) = fiber(i_node, 3) + (fiber_radius(i_node)) * ( n_prime(1,3) * cos(theta) + b_prime(1,3) * sin(theta));
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
          surfaceCorners, bundle_colors] = checkInputs(var,nvarsin)
%
%
% Preprocess the inputs tot he main function and return the fibers and all the default
% parameters.
%

% Check the number of inputs of the main function.
if (nvarsin < 1 || nvarsin > 8)
  error('[%s] Incorrect number of arguments (supplied: %i)\n',mfilename,nvarsin);
end
fiberRadius = 0.2;
% % % Reorganize the fibers.
% % %
% % % Strip off first and last points from strand
% % numFibers = length(var{1}.fibers);
% % fibers    = cell(numFibers,1);
% % parfor ll = 1:numFibers
% %   fibers{ll,1} = var{1}.fibers{ll}';
% %   %fibers{ll,2} = ll - 1;
% %   %fibers{ll,4} = ll - 1;
% %   %fibers{ll,3} = fiberRadius;
% % end
% % collection_name = var{1}.name;
% % 
% % % Figure handle passed in by fePlot
% % if (nvarsin < 2), fig_handle = figure;
% % else   fig_handle = var{2};
% % end

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
surfaceCorners = tube_divs + 1;

end % end function


