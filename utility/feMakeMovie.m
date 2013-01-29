function feMakeMovie(h,views,fileName)
% function feMakeMovie(algo,fibers,fileName)
%
% Make movies of various fiber gorups in a brain.
%
% Franco (c) Stanford Vista Team 2012

% ---- 
hold on
axis vis3d
set(h,'Position',[0 0 .25 .75],'Color',[0 0 0]);

vidObj = VideoWriter(fileName);
vidObj.FrameRate = 20;
vidObj.Quality   = 90;
open(vidObj);

% Rotate forward
for fi = 1:floor(views/2)
  camorbit(5,0,'data',[0 0 1])
  drawnow
  writeVideo(vidObj,getframe(h))
end
% Rotate back
for fi = 1:ceil(views/2)
  camorbit(-5,0,'data',[0 0 1])
  drawnow
  writeVideo(vidObj,getframe(h))
end

% Rotate back
for fi = 1:ceil(views/2)
  camorbit(-5,0,'data',[0 0 1])
  drawnow
  writeVideo(vidObj,getframe(h))
end
% Rotate forward
for fi = 1:floor(views/2)
  camorbit(5,0,'data',[0 0 1])
  drawnow
  writeVideo(vidObj,getframe(h))
end
fprintf('[%s] Writing movie...\n%s.avi\n',mfilename,fileName)
close(vidObj); 

