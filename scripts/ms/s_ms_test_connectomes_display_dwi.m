function s_ms_test_connectomes_display_dwi(trackingType,lmax,bval,rep)
%
% s_ms_test_connectomes_display_dwi(trackingType,lmax,bval,rep)
%
% Load a series of FE structures for conectomes of the occipital lobe.
% Saves figures of the occipital diffusion signal (measured, predicted, error).
%
% This is part of a series of reproducible science scritps to be published
% with the LiFE mansucript.
%
% Franco (c) Stanford Vista Team 2013

if notDefined('trackingType'),trackingType = 'deterministic';end
if notDefined('lmax'),        lmax         = 8;end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = 1;end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_figures');end
dirs = 100;
volume = [0 0 -3 dirs];

xlim = [0 50];
ylim = [-100 -60];

% High-resolution Anatomy
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';

% Information on the path to the files to load.
% This is where the inputs will be loaded from
feFileToLoad = msBuildFeFileName(trackingType,lmax,bval,rep,diffusionModelParams);

% Get the fe structure
load(feFileToLoad);

% Get the signal into an image.
dSig    = feGet(fe,'dsigdemeanedvox');
dSigImg = feValues2volume(dSig,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);

% Create the nifti structure
niM  = niftiCreate('data',dSigImg, ...
    'qto_xyz',feGet(fe,'xform img 2 acpc'), ...
    'fname','dwi_measured_signal', ...
    'data_type',class(dSigImg));

% Resample at the resolution of the T1.
%niM = mrAnatResampleToNifti(niM,t1File);

% Dispaly the good fibers 
fh = mrvNewGraphWin('Measured DW signal');
feDisplayBrainSlice( niM, volume);
view(0,90)
axis equal
hold on
set(gca,'xlim',[0 50],'ylim',[-100 -60])
set(fh,'Position',[0 0 .45 .95],'Color',[0 0 0]);
drawnow

saveFig(fh,fullfile(saveDir,sprintf('dwi_measured_signal_dir%i_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',dirs,trackingType,lmax,bval,rep, ...
  100*diffusionModelParams(1),100*diffusionModelParams(2))));

% Get the signal into an image.
dSig    = feGetRep(fe,'dsigdemeanedvox');
dSigImg = feValues2volume(dSig,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);

% Create the nifti structure
niM  = niftiCreate('data',dSigImg, ...
    'qto_xyz',feGet(fe,'xform img 2 acpc'), ...
    'fname','dwi_measured_signal', ...
    'data_type',class(dSigImg));

% Resample at the resolution of the T1.
%niM = mrAnatResampleToNifti(niM,t1File);

% Dispaly the good fibers 
fh = mrvNewGraphWin('Measured DW signal 2');
feDisplayBrainSlice( niM, volume);
view(0,90)
axis equal
hold on
set(gca,'xlim',[0 50],'ylim',[-100 -60])
set(fh,'Position',[0 0 .45 .95],'Color',[0 0 0]);
drawnow

saveFig(fh,fullfile(saveDir,sprintf('dwi_measured_signal_2_dir%i_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',dirs,trackingType,lmax,bval,rep, ...
  100*diffusionModelParams(1),100*diffusionModelParams(2))));

%% Predicted signal
% Get the signal into an image.
pSig    = feGet(fe,'psigfvox');
pSigImg = feValues2volume(pSig,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);
keyboard
% Create the nifti structure
niP  = niftiCreate('data',pSigImg, ...
    'qto_xyz',feGet(fe,'xform img 2 acpc'), ...
    'fname','dwi_predicted_signal', ...
    'data_type',class(pSigImg));

% Resample at the resolution of the T1.
%niP = mrAnatResampleToNifti(niP,t1File);

% Dispaly the good fibers 
fh = mrvNewGraphWin('Predicted DW signal');
feDisplayBrainSlice( niP, volume);
view(0,90)
axis equal
hold on
set(gca,'xlim',[0 50],'ylim',[-100 -60])
set(fh,'Position',[0 0 .45 .95],'Color',[0 0 0]);
drawnow

saveFig(fh,fullfile(saveDir,sprintf('dwi_predicted_signal_dir%i_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',dirs,trackingType,lmax,bval,rep, ...
  100*diffusionModelParams(1),100*diffusionModelParams(2))));

%% Error signal
% Get the signal into an image.
eSig    = feGet(fe,'voxrmse');%'resfibervox')';
eSigImg = feValues2volume(eSig,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);

% Create the nifti structure
ni  = niftiCreate('data',eSigImg, ...
    'qto_xyz',feGet(fe,'xform img 2 acpc'), ...
    'fname','dwi_predicted_signal', ...
    'data_type',class(eSigImg));

% Resample at the resolution of the T1.
%ni = mrAnatResampleToNifti(ni,t1File);

% Dispaly the good fibers 
fh = mrvNewGraphWin('Predicted DW signal');
feDisplayBrainSlice( ni, volume);
view(0,90)
axis equal
hold on
set(gca,'xlim',[0 50],'ylim',[-100 -60])
set(fh,'Position',[0 0 .45 .95],'Color',[0 0 0]);
drawnow

saveFig(fh,fullfile(saveDir,sprintf('dwi_rmse_signal_dir%i_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',dirs,trackingType,lmax,bval,rep, ...
  100*diffusionModelParams(1),100*diffusionModelParams(2))));

end

%---------------------------------------%
function makeAllMapsLocal(fe,tag,mapsdir,movieFileName,  dwiNifti)

xform = dwiNifti.qto_xyz;

% Measured
% Get the diffusion signal
mapName = fullfile(mapsdir,  [movieFileName, '_arcuate_dw_signal',tag]);

dSig    = feGet(fe,'voxdsigfull')';
dSigImg = feValues2volume(dSig,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);


% Predicted
mapName = fullfile(mapsdir,  [movieFileName, '_predicted_dw_signal',tag]);
pSig    = feGet(fe,'psigfvoxelwisebyvoxel');
pSigImg = feValues2volume(pSig,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);


% GEt the residual signal
mapName = fullfile(mapsdir,  [movieFileName, '_residual_dw_signal',tag]);
res     = feGetRep(fe,'voxressigfullvoxelwise');
resImg  = feValues2volume(res,feGet(fe,'roicoords'),feGet(fe,'volumesize')-[0 0 0 10]);

% Get the rRMSE
rr      = feGetRep(fe,'voxrmseratiovoxelwise');
rrImg   = feValues2volume(rr,feGet(fe,'roicoords'),feGet(fe,'mapsize'));
mapName = fullfile(mapsdir,['arcuate_rmse_ratio',tag]);

% Get the RMSE
rmse    = feGetRep(fe,'voxrmsevoxelwise');
rmseImg = feValues2volume(rmse,feGet(fe,'roicoords'),feGet(fe,'mapsize'));
mapName = fullfile(mapsdir,['arcuate_rmse',tag]);
%feWriteValues2nifti(rmseImg,mapName,xform);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(h),figName));

end