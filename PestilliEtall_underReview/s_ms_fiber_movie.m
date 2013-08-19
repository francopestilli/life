function s_ms_fiber_movie(savemovie,moviedir)
%
% Script used to visualize the tensor as a 3D and 3D signal.
% Saves movies
%
% Franco
global gPlotBox
if notDefined('moviedir'),moviedir   = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/fiber_ellipsoid_movies';end
if notDefined('savemovie'),savemovie = 1;end
gPlotBox = 0;

% PARAMETERS
diffusionModelParams = [1,0];       % The parameters of the tensor model AD, RD
sdCutoff             = [3.3];     % One per conenctome, generally smaller for smaller lmax values
nSkip = 6; % How many fiber nodes to skip at each step of the tensor
fibersIdx = {2200,2300,2330};
fcolor = {[.95 .7 .8],[.95 .4 .5],[.92 .1 .15]};
nframes = 150;
nfibers = 320;

% DIRECTORY TO LOAD FILES FROM:
% DWI data
dataRootPath  = fullfile('/biac2','wandell6','data','frk','life_dti','FP20120420');
subfolders    = fullfile('150dirs_b1000_1');
baseDir       = fullfile(dataRootPath,subfolders);
dtFile        = fullfile(baseDir,'dt6.mat');
dwiFile       = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
dwiFileRepeat = fullfile(dataRootPath,'raw','0011_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
t1File        = fullfile(dataRootPath,'t1','t1.nii.gz');

% ROIs connectomes
projectDir           = '/azure/scr1/frk/150dirs_b1000_b2000_b4000';
connectSubfolders    = {'life_mrtrix_rep1'};%,'life_mrtrix_rep2','life_mrtrix_rep3'};
connectomeFile       = { '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax12_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_prob-500000.pdb'};

% DIRECTORIES and FILES TO SAVE:
saveDir              =   fullfile(projectDir,'results');
fe_structuresDir     =  'fe_arcuate_importance';
movieFileName        = 'arcuate_fiber_ellipsoid';

irep = 1;
ii = 1;

% Set up file names
% Find an identificative name for the connectome that is shortenough:
cName         = [connectomeFile{ii}(1:57),'_',connectomeFile{ii}(end-16:end-4)];

% Name and path for savign the fe structures
feSaveDir     = fullfile(saveDir,connectSubfolders{irep},fe_structuresDir);

% Build the full name of the two fascicles FE's structures
feSaveNameAll     = sprintf('%s_diffModAx%sRd%s_%s',cName,num2str(100*diffusionModelParams(1)), ...
  num2str(100*diffusionModelParams(2)));
feSaveNameArcuate = sprintf('arcuateFE_%s_sd%2.0f',feSaveNameAll,100*sdCutoff(1));

%% Build the LiFE model around the volume of the Left Arcuate Fasciculum
if ~(exist(fullfile(feSaveDir,[feSaveNameArcuate,feSaveNameArcuate,'.mat']),'file') == 2)
  fprintf('[%s] CANNOT FIND the LiFE model for the arcuate only \n%s\n ======================================== \n\n',mfilename,feSaveNameArcuate)
  keyboard
else
  fprintf('[%s] FOUND fe FILE ofr acruate ONLY (no connectome), Loading it: \n%s\n ======================================== \n\n',mfilename,feSaveNameArcuate)
  load(fullfile(feSaveDir,[feSaveNameArcuate,feSaveNameArcuate,'.mat']),'fe');
end
disp('Loading T1')
nifti = readFileNifti('/biac2/wandell6/data/arokem/ModelFits/FP20120420/t1/t1.nii.gz');

disp('MAKING MOVIE...')

% Open a video object
if savemovie
  vidObj = VideoWriter(fullfile(moviedir,  movieFileName));
  vidObj.FrameRate = 10;
  vidObj.Quality = 90;
  open(vidObj);
end

% Get the fibers
fg = feGet(fe,'fibers acpc'); clear fe
h = figure;

fgW=fg;
% Display a large numer of fibers then smaller and smaller up to 1
fiber_indices = randsample(1:length(fg.fibers),nfibers);
fibers_to_show = [length(fiber_indices) ceil(length(fiber_indices)/2) ceil(length(fiber_indices)/4) ...
                    ceil(length(fiber_indices)/8) ceil(length(fiber_indices)/16) ceil(length(fiber_indices)/32) ...
                    ceil(length(fiber_indices)/64)];
for fib = fibers_to_show
  fgW.fibers = fg.fibers(fiber_indices(1:fib));
  drawnow
  [h,sh] = feConnectomeDisplay(fgW,h,[.9 .7 .7]);
  set(gcf,'Position',[0 0 0.55 0.95]);
  if fib ==fibers_to_show(1)
    axis off
    set(gcf,'Color',[0 0 0]);
    axis tight
    view(-90,0);
    set(gca,'looseInset',get(gca,'TightInset'))
    ch = camlight('left');
    hold on
    mctDisplayBrainSlice(nifti,[-15 0 0]);
    
    % Save a figure
    feSavefig(gcf,'figName','arcuate_fiber_movie_fig.jpg','figDir',moviedir,'figType','jpg')
    feSavefig(gcf,'figName','arcuate_fiber_movie_fig.png','figDir',moviedir,'figType','png')
  end
  drawnow

  % Make a short segment wtih the first fibe ronly
  if savemovie
    for jj = 1:10 % 5 frames per step
      writeVideo(vidObj,getframe(h))
    end
  end
  delete(sh);
end
delete(ch)

fg1 =fg;
fg1.fibers = {};

for ii =1:length(fibersIdx)
  fg1.fibers{ii} = fg.fibers{fibersIdx{ii}};
end
thisFiber = fg1.fibers{1};
fgtmp = fg1;fgtmp.fibers = {};
fgtmp.fibers{1} = thisFiber;

% Calculate a simulated tensor for each point and each fiber.
% The tensors all have a common axial and radial diffusivity
d_ad = 3; d_rd = 0.8;
dParms(1) = d_ad; dParms(2) = d_rd; dParms(3) = d_rd;
fg1.tensors = feComputeCanonicalDiffusion(fg1.fibers,dParms);

% plot  fiber
feConnectomeDisplay(fgtmp,h,fcolor{1});
set(gcf,'Position',[0 0 0.55 0.95]);
ch = camlight('left');
hold on
drawnow

fgtmp = fg1;
fgtmp.fibers = {};
for iff = 1:length(fibersIdx)
  thisFiber = fg1.fibers{iff};
  fgtmp = fg1;fgtmp.fibers = {};
  fgtmp.fibers{1} = thisFiber;
  thisTensor  = fg1.tensors{iff};
  
  if iff >1
    [h,sh] = feConnectomeDisplay(fgtmp,h,fcolor{iff});  
    set(gcf,'Position',[0 0 0.55 0.95]);
    %ch(iff-1) = camlight('left');
    hold on
    drawnow
  end
  
  for ii=[size(thisFiber,2):-nSkip:1]
    Q = reshape(thisTensor(ii,:),3,3);
    C = thisFiber(:,ii)';
    hs = ellipsoidFromTensorLocal(Q,C,5);
    drawnow
    if savemovie
      for jj = 1:5 % 5 frames per step
      writeVideo(vidObj,getframe(gcf))
      end
    end
    delete(hs);    
  end
end

% Make a short segment wtih the three fiber 
if savemovie
  for ii=nframes
    writeVideo(vidObj,getframe(gcf))
  end
end
delete(ch);

fgW=fg;
% Display a large numer of fibers then smaller and smaller up to 1
for fib = fliplr(fibers_to_show)
  fgW.fibers = fg.fibers(fiber_indices(1:fib));
 
  [h,sh] = feConnectomeDisplay(fgW,h,[.9 .7 .7]);
  set(gcf,'Position',[0 0 0.55 0.95]);
  ch = camlight('left');
  hold on
  drawnow

  % Make a short segment wtih the first fibe ronly
  if savemovie
    for jj = 1:10 % 5 frames per step
      writeVideo(vidObj,getframe(h))
    end
  end
  if fib < 320
    delete(sh);
  end
  delete(ch);
end
    
if savemovie
  close(vidObj);
end

disp('DONE MAKING MOVIE...')

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
if notDefined('N'), N = 120; end      % number of points on the sphere
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

hs = surf(x,y,z,repmat(256,size(z)),'EdgeAlpha',0.001,'FaceAlpha',0.85);

cmap = hot(255);
colormap([cmap; .25 .35 .65])
end
