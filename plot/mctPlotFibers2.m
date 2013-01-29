function fh = mctPlotFibers(fg,figName,numClones,mmPerNode)

% Make a 3d plot fo the fibers
c = {'r','g','b','c','m','k','y'};
scale = 0.0125; % plot the tensor at 10% the original size
if mmPerNode > 1
    nSkip = 1;
else
   nSkip = 4;
end

newFigure = 0; % plot ont he current fibure

main_fig = mrvNewGraphWin(figName);
set(main_fig,'Units','normalized')
set(main_fig, 'Position', [0.3 0.25 0.4 0.5]);
set(main_fig, 'DoubleBuffer', 'on');
cameratoolbar('Show');
cameratoolbar('SetMode','orbit');
clf;
whitebg(main_fig,'black');
hold on;
    
for ii =1:length(fg.fibers)
    % use mod or rem to switch colors for fiber groups    if mod(numFibers,numClones)
    hold on
    fh = plot3(fg.fibers{ii}(1,:),fg.fibers{ii}(2,:),fg.fibers{ii}(3,:),'.-','color','b','lineWidth',2);
    
    % Draw Vectors directions (quivers) of x and x + gradient(x).
    % This shows the points and the tangent to the curve
    % The final 0 is essential for handling Matlab's little arrowheads.
    imgGradient = gradient(fg.fibers{ii});
    
    quiver3(fg.fibers{ii}(1,:),fg.fibers{ii}(2,:),fg.fibers{ii}(3,:),...
            imgGradient(1,:),imgGradient(2,:),imgGradient(3,:),0,'r-')
 
    if isfield(fg,'Q')
        for in=1:nSkip:size(fg.fibers{ii},2)
            Q = reshape(fg.Q{ii}(in,:),3,3) * scale;
            C = fg.fibers{ii}(:,in)' ;
            ellipsoidFromTensor(Q,C,31,newFigure);
            hold on;
        end
    end
end

axis([.5 5.5 .5 5.5 .5 1.5])
set(gca,'yTick', .5 + [1 2 3 4 5], 'xTick', .5 + [1 2 3 4 5])
grid on
hold off;

% From display_strands
set(gcf, 'color', [0 0 0]);
% 	xlabel('X-axis', 'color', [1 1 1]);
% 	ylabel('Y-axis', 'color', [1 1 1]);
% 	zlabel('Z-axis', 'color', [1 1 1]);

label_handle = get(gca,'xlabel');
set(label_handle, 'string', 'X-axis', 'color', [1 1 1]);
% 	label_position = get(label_handle, 'position');
% 	set(label_handle, 'position', (label_position + [0,  0, 0]));
%
label_handle = get(gca,'ylabel');
set(label_handle, 'string', 'Y-axis', 'color', [1 1 1]);
% 	label_position = get(label_handle, 'position');
% 	set(label_handle, 'position', (label_position + [0,  0, 0]));

label_handle = get(gca,'zlabel');
set(label_handle, 'string', 'Z-axis', 'color', [1 1 1]);
% 	label_position = get(label_handle, 'position');
% 	set(label_handle, 'position', (label_position + [0,  0, 0]));




h = get (gca, 'children');
daspect ([ 1 1 1 ]);
if (plot_3D)
    set(h,'edgecolor', 'none')
end
light
lighting gouraud;


function [t,n,b]=frame(x,y,z,vec)

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
end