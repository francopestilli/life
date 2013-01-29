function hRoi = mctDisplayRoi(fh,roi)
% function hs = mctDisplayRoi(roi)
%
% Plots an roi in as a semi-transparent surface in 3D matlab plot.
%
% (c) Franco Stanford VISTA Team

rendering = 'tri';
figure(fh);

switch rendering
  case 'tetra'
    T = DelaunayTri(roi.coords(:,1),roi.coords(:,2),roi.coords(:,3));
    hRoi = tetramesh(T,'Marker','.', 'MarkerFaceColor',[.8 .6 .2], ...
                       'FaceAlpha',0.35,'FaceLighting','gouraud', ...
                       'EdgeAlpha',0.001,'EdgeColor',[.5 .5 .5],'EdgeLighting','phong', ...
                       'SpecularColorReflectance',.5, 'MarkerSize',1);
    
  case 'tri'
    [~,S] = alphavol(roi.coords,1.8);
    hRoi  = trisurf(S.bnd, roi.coords(:,1),roi.coords(:,2),roi.coords(:,3),'EdgeAlpha',0,'FaceAlpha',1,'LineWidth',0.00001);
    set(hRoi, 'FaceColor', [.9 .25 .1])
    light, lighting phong;
  otherwise
    keyboard
end
