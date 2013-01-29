function fh = mctPlotFibers(fg,figName,whichFibers,mmPerNode)

% Make a 3d plot fo the fibers
c = {'r','g','b','c','m','k','y'};
scale = 0.05; % plot the tensor at scale*100 = % the original size
if mmPerNode > 1
  nSkip = 1;
else
  nSkip = ceil(mmPerNode * .125);
end

newFigure = 0; % plot ont he current fibure
mrvNewGraphWin(figName);

for ii =1:length(fg.fibers)
  hold on
  fh = plot3(fg.fibers{ii}(1,:),fg.fibers{ii}(2,:),fg.fibers{ii}(3,:),'o-','color','b','lineWidth',2);
  
  % Draw Vectors directions (quivers) of x and x + gradient(x).
  % This shows the points and the tangent to the curve
  % The final 0 is essential for handling Matlab's little arrowheads.
  %imgGradient = gradient(fg.fibers{ii});
  
  %quiver3(fg.fibers{ii}(1,:),fg.fibers{ii}(2,:),fg.fibers{ii}(3,:),...
  %        imgGradient(1,:),imgGradient(2,:),imgGradient(3,:),0,'r-')
  
  if isfield(fg,'Q')
    for in = 1 : nSkip : size(fg.fibers{ii},2)
      Q = reshape(fg.Q{ii}(in,:),3,3) * scale;
      C = fg.fibers{ii}(:,in)' ;
      ellipsoidFromTensor(Q,C,20,newFigure);
      hold on;
    end
  end
end

axis([.5 5.5 .5 5.5 .5 1.5])
set(gca,'yTick', .5 + [1 2 3 4 5], 'xTick', .5 + [1 2 3 4 5])
grid on
axis equal


