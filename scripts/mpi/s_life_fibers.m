function s_life_fibers(savemovie,moviedir,showRoi)
%
% Script used to visualize the tensor as a 3D and 3D signal.
% Saves movies
%
% Franco
global gPlotBox
if notDefined('algo');   algo = 'TEND';end % TEND or STT
if notDefined('this_roi');   this_roi = 'mct_roi_RTO_12.mat';end% 'mct_r_lgn.mat';end
if notDefined('moviedir'),moviedir   = '/biac2/wandell6/data/frk/LiFE/results/mpi/fiber_movies';end
if notDefined('savemovie'),savemovie = 1;end
if notDefined('showRoi'), showRoi = 0; end
gPlotBox = 0;

dataRootPath = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
subfolders   = fullfile('150dirs_b1000_1');
dataDir      = fullfile(dataRootPath,subfolders);
dtFile       = fullfile(dataDir,'dt6.mat');


roiDir   = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/To12_test';
roiName  = fullfile(roiDir,this_roi);
roi      = dtiReadRoi(roiName);

fiberDir     = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/To12_test';
fgName       = fullfile(fiberDir,sprintf('mct_fg_RTO_12_%s.mat',algo));
fg     = dtiLoadFiberGroup(fgName);

disp('Loading T1')
nifti = readFileNifti('/biac2/wandell6/data/arokem/ModelFits/FP20120420/t1/t1.nii.gz');

fibers = fg.fibers;
fg.fibers = {};
fg.fibers{1} = fibers{3};

fh = mctNfgDisplayStrands(fg);
keyboard
hold on
hm = mctDisplayBrainSlice(nifti,[0 0 -7]);
if showRoi, mctDisplayRoi(fh,roi); end
keyboard
set(gcf,'Position',[0 0 0.9 0.9]);
axis off
set(gcf,'Color',[0 0 0]);
view(161,30)
set(gca,'looseInset',get(gca,'TightInset'))
axis tight


% For each point on each fiber, we make its tensor.
thisFiber = fg.fibers{1};

% Calculate a simulated tensor for each point and each fiber.
% The tensors all have a common axial and radial diffusivity
d_ad = 1.7; d_rd = 0.2;
dParms(1) = d_ad; dParms(2) = d_rd; dParms(3) = d_rd;
fg.tensors = fgTensors(fg,dParms);  % Should become fgSet
thisTensor    = fg.tensors{1};


% Open a video object
if savemovie
  vidObj = VideoWriter(fullfile(moviedir,  sprintf('%s','mct_fiber_movie_ilf_no_roi_%s',algo)));
  vidObj.FrameRate = 7.5;
  vidObj.Quality = 100;
  open(vidObj);
end


nSkip = 5;
c = 1;  
for ii=[1:nSkip:size(thisFiber,2), size(thisFiber,2):-nSkip:1]
    Q = reshape(thisTensor(ii,:),3,3);
    C = thisFiber(:,ii)';
    hs = ellipsoidFromTensorLocal(Q,C,5); 
    drawnow
    
    if savemovie 
      m(c) = getframe(gcf);      
      writeVideo(vidObj,m(c))
    end
    delete(hs);
    c = c + 1;
end

if savemovie
  close(vidObj);
end

% This is an example of how to compute the AD, RD and PDD 
% Calculate eigenvalues and eigenvectors
%[V,D] = eigs(Q);
%val = diag(D);
%ad = val(1);
%rd = (val(2) + val(3))/2;
%pdd = V(:,1);
%ad, rd
%pdd

end

%%%%%%%%%%%%%
% makeMovie %
%%%%%%%%%%%%%
function m = makeMovie(fh,ah,filename,replay)
if notDefined('replay'),replay=0;end

% rotate the sphere and collect frames.
c = 1;
rotations = [10 12 14 16 20 25 30:10:70 75 80 84 86 88 90];
rotations = [rotations, fliplr(rotations)];
rotations = (-13:.4:13).^2 + 1;
%rotations = [logspace(-3,0,20) * 180, fliplr(logspace(-3,0,20)) * 180];
set(fh,'Position',[0 0 1600 1200],'Visible','on');
axis off;

for ii = rotations 
   view(ah,ii,10);
   drawnow
   m(c) = getframe(fh);
   c = c + 1;
end


if replay
  % use 1st frame to get dimensions
  mfh       = figure;
  [h, w, p] = size(m(end).cdata);
  set(mfh, 'Position', [10 10 w h]); 
  movie(mfh, m,1,12,[0 0 0 0]);
end

% Save movie to disk
movie2avi(m, sprintf('%s.avi',filename),'Compression','None');
end


%%%%%%%%%%%%%%%%%%%%%%%
% ellipsoidFromTensor %
%%%%%%%%%%%%%%%%%%%%%%%
function hs = ellipsoidFromTensorLocal(Q,C,tensorDim)
%Plot a diffusion ellipsoid from a tensor
%
% This is a copy of ellipsoid for tensor for local edits.
%
% (c) Stanford VISTA Team
if notDefined('N'), N = 60; end      % number of points on the sphere
if notDefined('C'), C = [0,0,0]; end % coordinates of the points in 3D space

% Make unit vectors on a sphere, u
[x,y,z] = sphere(N);   % build a sphere sampled with N resolution
u = [x(:), y(:), z(:)];% reorganize the univectors on the sheres (basically the diffusion directions scaled to the unit sphere)

sphAdc = zeros(size(u,1),1);
for ii=1:size(u,1)
  % Notice that this is multiplying by inv(Q), not Q.
  % Uncomment to see the equivalence of
  % sphADC(ii) = u(ii,:)*inv(Q)*u(ii,:)';
  sphAdc(ii) = u(ii,:)*(Q\u(ii,:)');
end

% The sqrt ADC is multipled into the unit vector because the factor is
% multiplied in twice (left and right side of inv(Q)).
v = diag(1 ./ sphAdc.^0.5)*u;  %Unit vectors in the rows.  Each row scaled.

% Now reshape and plot
sz = size(x);
x = reshape(v(:,1),sz) .* tensorDim;
y = reshape(v(:,2),sz) .* tensorDim;
z = reshape(v(:,3),sz) .* tensorDim;
x = x + C(1); y = y + C(2); z = z + C(3);

hs = surf(x,y,z,repmat(256,size(z)),'EdgeAlpha',0.01);
cmap = hot(255);
colormap([cmap; .25 .35 .65])
end
