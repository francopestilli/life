function s_ms_fiber_dwi_prediction_movie(slice,savemovie,moviedir)
%
% Script used to visualizepredicted signal from a fascicle.
% Saves movies
%
% Franco
global gPlotBox
if notDefined('moviedir'), moviedir  = '~/Dropbox/arcuate_cst_movies/signal';end
if notDefined('mapsdir'),  mapsdir   = '~/Dropbox/arcuate_cst_maps/';end
if notDefined('savemovie'),savemovie = 1;end
gPlotBox = 0;

% PARAMETERS
diffusionModelParams = [1,0];       % The parameters of the tensor model AD, RD
sdCutoff             = [3.2];     % One per conenctome, generally smaller for smaller lmax values
nframes = 12;

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
saveDir              = fullfile(projectDir,'results');
fe_structuresDir     = 'fe_arcuate_cst_test';
movieFileName        = 'arcuate_fiber_ellipsoid';

irep = 1;ii = 1;

% Set up file names
% Find an identificative name for the connectome that is shortenough:
cName         = [connectomeFile{ii}(1:57),'_',connectomeFile{ii}(end-16:end-4)];

% Name and path for savign the fe structures
feSaveDir     = fullfile(saveDir,connectSubfolders{irep},fe_structuresDir);

% Build the full name of the two fascicles FE's structures
feSaveNameAll     = sprintf('%s_diffModAx%sRd%s_%s',cName,num2str(100*diffusionModelParams(1)), ...
  num2str(100*diffusionModelParams(2)));
feSaveNameArcuate = 'arcuateFE_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax12__m_prob-500000_diffModAx100Rd0_arcuateFE_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax12__m_prob-500000_diffModAx100Rd0_SD320';
sprintf('arcuateFE_%s_sd%2.0f',feSaveNameAll,100*sdCutoff(1));

%% Build the LiFE model around the volume of the Left Arcuate Fasciculum
fprintf('[%s] LOAD fe FILE for arcuate ONLY: \n%s\n ======================================== \n\n',mfilename,feSaveNameArcuate)
load(fullfile(feSaveDir,[feSaveNameArcuate,'.mat']),'fe');

disp('Loading DWI')
dwiNifti = readFileNifti(dwiFile);

% Make a movie
makeAllMapsLocal(fe,'_arcuate_only',mapsdir,movieFileName,  dwiNifti)

%% Now load the Arcuate PLUS CST strucutre and re do movies and figures
feSaveNameArcuate = 'arcuateCstUnionFE_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax12__m_prob-500000_diffModAx100Rd0_arcuateCstUnionFE_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax12__m_prob-500000_diffModAx100Rd0_SD320';

%% Build the LiFE model around the volume of the Left Arcuate Fasciculum
fprintf('[%s] LOAD fe FILE for arcuate AND CST: \n%s\n ======================================== \n\n',mfilename,feSaveNameArcuate)
load(fullfile(feSaveDir,[feSaveNameArcuate,'.mat']),'fe');

% Make a movie
makeAllMapsLocal(fe,'_arcuate_cst',mapsdir,movieFileName,  dwiNifti)

end

%---------------------------------------%
function makeAllMapsLocal(fe,tag,mapsdir,movieFileName,  dwiNifti)
disp('making movies...')

xform = dwiNifti.qto_xyz;

% Measured
if strcmpi(tag,'_arcuate_cst')    
  % Get the diffusion signal
  mapName = fullfile(mapsdir,  [movieFileName, '_arcuate_dw_signal',tag]);
  
  dSig    = feGet(fe,'voxdsigfull')';
  dSigImg = feValues2volume(dSig,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);
  feWriteValues2nifti(dSigImg,mapName,xform);
end

% Predicted
mapName = fullfile(mapsdir,  [movieFileName, '_predicted_dw_signal',tag]);
pSig    = feGet(fe,'psigfvoxelwisebyvoxel');
pSigImg = feValues2volume(pSig,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);
feWriteValues2nifti(pSigImg,mapName,xform);


% GEt the residual signal
mapName = fullfile(mapsdir,  [movieFileName, '_residual_dw_signal',tag]);
res     = feGetRep(fe,'voxressigfullvoxelwise');
resImg  = feValues2volume(res,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);
feWriteValues2nifti(resImg,mapName,xform);

% Get the rRMSE
rr      = feGetRep(fe,'voxrmseratiovoxelwise');
rrImg   = feValues2volume(rr,feGet(fe,'roicoords'),feGet(fe,'mapsize'));
mapName = fullfile(mapsdir,['arcuate_rmse_ratio',tag]);
feWriteValues2nifti(rrImg,mapName,xform);

% Get the RMSE
rmse    = feGetRep(fe,'voxrmsevoxelwise');
rmseImg = feValues2volume(rmse,feGet(fe,'roicoords'),feGet(fe,'mapsize'));
mapName = fullfile(mapsdir,['arcuate_rmse',tag]);
feWriteValues2nifti(rmseImg,mapName,xform);

end


%---------------------------------------%
function makeAllMoviesLocal(fe,tag,moviedir,movieFileName, slice, nframes,savemovie,  dwiNifti)
disp('making movies...')

if strcmpi(tag,'_arcuate_cst')
  % Show a slice
  h = figure('name','Original DW Signal','color','k');
  set(gcf,'Units','normalized','Position',[0 0 0.2 0.7]);
  % Open a video object
  if savemovie
    vidObj = VideoWriter(fullfile(moviedir,  [movieFileName, '_original_dw_signal_slice',num2str(slice),tag]));
    vidObj.FrameRate = 10;
    vidObj.Quality = 90;
    open(vidObj);
  end
  
  img = flipud(squeeze(dwiNifti.data(:,:,slice,11))');
  imagesc(img);drawnow
  axis equal
  axis off
  colormap bone
  hold on
  for ii =1:5:150-1
    img = flipud(squeeze(dwiNifti.data(:,:,slice,11+ii))');
    imagesc(img);drawnow
    % Make a short segment wtih the first fibe ronly
    if savemovie
      for jj = 1:nframes % 5 frames per step
        writeVideo(vidObj,getframe(gcf))
      end
    end
  end
  if savemovie
    close(vidObj);
  end
  
  
  % Get the diffusion signal
  h = figure('name','Original DW Signal Arcuate','color','k');
  set(gcf,'Units','normalized','Position',[0 0 0.2 0.7]);
  % Open a video object
  if savemovie
    vidObj = VideoWriter(fullfile(moviedir,  [movieFileName, '_arcuate_dw_signal_slice',num2str(slice),tag]));
    vidObj.FrameRate = 10;
    vidObj.Quality = 90;
    open(vidObj);
  end
  
  dSig = feGet(fe,'voxdsigfull')';
  dSigImg = feValues2volume(dSig,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);
  
  img = flipud(squeeze(dSigImg(:,:,slice,1))');
  imagesc(img);drawnow
  axis equal
  axis off
  colormap bone
  hold on
  for ii =1:5:150-1
    img = flipud(squeeze(dSigImg(:,:,slice,1+ii))');
    imagesc(img);drawnow
    % Make a short segment wtih the first fibe ronly
    if savemovie
      for jj = 1:nframes % 5 frames per step
        writeVideo(vidObj,getframe(gcf))
      end
    end
  end
  if savemovie
    close(vidObj);
  end
end

% Get the predicted signal
h = figure('name','Predicted DW Signal','color','k');
  set(gcf,'Units','normalized','Position',[0 0 0.2 0.7]);
% Open a video object
if savemovie
  vidObj = VideoWriter(fullfile(moviedir,  [movieFileName, '_predicted_dw_signal_slice',num2str(slice),tag]));
  vidObj.FrameRate = 10;
  vidObj.Quality = 90;
  open(vidObj);
end

pSig = feGet(fe,'psigfvoxelwisebyvoxel');
pSigImg = feValues2volume(pSig,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);
img = flipud(squeeze(pSigImg(:,:,slice,1))');
imagesc(img);drawnow
axis equal
axis off
colormap copper
hold on
for ii =1:5:150-1
  img = flipud(squeeze(pSigImg(:,:,slice,1+ii))');
  imagesc(img);drawnow
  % Make a short segment wtih the first fibe ronly
  if savemovie
    for jj = 1:nframes % 5 frames per step
      writeVideo(vidObj,getframe(gcf))
    end
  end
end
if savemovie
  close(vidObj);
end

% GEt the residual signal
h = figure('name','Residual DW Signal','color','k');
  set(gcf,'Units','normalized','Position',[0 0 0.2 0.7]);
% Open a video object
if savemovie
  vidObj = VideoWriter(fullfile(moviedir,  [movieFileName, '_residual_dw_signal_slice',num2str(slice),tag]));
  vidObj.FrameRate = 10;
  vidObj.Quality = 90;
  open(vidObj);
end

res  = feGetRep(fe,'voxressigfullvoxelwise');
resImg = feValues2volume(res,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);
img = flipud(squeeze(resImg(:,:,slice,1))');
imagesc(img);drawnow
axis equal
axis off
colormap hot
hold on
for ii =1:5:150-1
  img = flipud(squeeze(resImg(:,:,slice,1+ii))');
  imagesc(img);drawnow
  % Make a short segment wtih the first fibe ronly
  if savemovie
    for jj = 1:nframes % 5 frames per step
      writeVideo(vidObj,getframe(gcf))
    end
  end
end
if savemovie
  close(vidObj);
end

% Get the rRMSE
rr  = feGetRep(fe,'voxrmseratiovoxelwise');
rrImg = feValues2volume(rr,feGet(fe,'roicoords'),feGet(fe,'mapsize'));
h = figure('name','RMSE Ratio DW Signal','color','k');
  set(gcf,'Units','normalized','Position',[0 0 0.2 0.7]);
img = flipud(squeeze(rrImg(:,:,slice,1))');
imagesc(img);drawnow
axis equal
axis off
colormap hot
colorbar
% Save a figure
feSavefig(gcf,'figName',['arcuate_rmse_ratio_slice',num2str(slice),tag,'.jpg'],'figDir',moviedir,'figType','jpg')
feSavefig(gcf,'figName',['arcuate_rmse_ratio_slice',num2str(slice),tag,'.png'],'figDir',moviedir,'figType','png')
feSavefig(gcf,'figName',['arcuate_rmse_ratio_slice',num2str(slice),tag,'.eps'],'figDir',moviedir,'figType','eps')

% Get the RMSE
rmse  = feGetRep(fe,'voxrmsevoxelwise');
rmseImg = feValues2volume(rmse,feGet(fe,'roicoords'),feGet(fe,'mapsize'));
h = figure('name','RMSE  DW Signal','color','k');
  set(gcf,'Units','normalized','Position',[0 0 0.2 0.7]);
img = flipud(squeeze(rmseImg(:,:,slice,1))');
imagesc(img);drawnow
axis equal
axis off
colormap hot
colorbar 
% Save a figure
feSavefig(gcf,'figName',['arcuate_rmse_slice',num2str(slice),tag,'.jpg'],'figDir',moviedir,'figType','jpg')
feSavefig(gcf,'figName',['arcuate_rmse_slice',num2str(slice),tag,'.png'],'figDir',moviedir,'figType','png')
feSavefig(gcf,'figName',['arcuate_rmse_slice',num2str(slice),tag,'.eps'],'figDir',moviedir,'figType','eps')

end

