function s_ms_test_connectomes_display_fg(trackingType,lmax,bval,rep)
%
% s_ms_test_connectomes_display_fg(trackingType,lmax,bval,rep)
%
% Load a series of FE structures for conectomes of the occipital lobe.
% Saves figures of the occipital conectomes.
%
% This is part of a series of reproducible science scritps to be published
% with the LiFE mansucript.
%
% Franco (c) Stanford Vista Team 2013


if notDefined('trackingType'),trackingType = 'deterministic';end
if notDefined('lmax'),        lmax         = 16;end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = 1;end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_figures');end

% Perecentile cut-off for the weights
weightsPercentile = [80];
color = [.2 .6 .75];

% High-resolution Anatomy
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';

% Information on the path to the files to load.
% This is where the inputs will be loaded from
feFileToLoad = msBuildFeFileName(trackingType,lmax,bval,rep,diffusionModelParams);

% Get the fe structure
load(feFileToLoad);

% We will now show and save out three connectomes. The full one, the one
% with the best fibers and the one with only the fibers rejected by the
% life fit.

% Find the fibers that have non-zero weight for life (Hits)
% Rough way to the the voxel-wise weights. Think more what to do wit this.
w.all = sum(fe.life.voxfit.weights);

% Find the positive non-zero weights
w.good  = w.all > 0 ;

% Extract the good fibers
fg = fgExtract(feGet(fe,'fibers acpc'),find(w.good),'keep');

% Dispaly the good fibers 
fh = figure('name','Accepted fibers');
feConnectomeDisplay(feSplitLoopFibers( fg ),fh,[.3 .85 .4], [], [], .35)
hold on
mctDisplayBrainSlice( niftiRead(t1File), [0 0 -10]);
set(gca,'ylim',[-126 -55],'xlim',[-5 80])
set(fh,'Position',[0 0 .45 .95],'Color',[0 0 0]);
drawnow

saveFig(fh,fullfile(saveDir,sprintf('connectome_HIT_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',trackingType,lmax,bval,rep, ...
  100*diffusionModelParams(1),100*diffusionModelParams(2))));

% Extract the bad fibers
fg = fgExtract(feGet(fe,'fibers acpc'),find(~w.good),'keep');

% Dispaly the bad fibers 
fh = figure('name','Rejected fibers');
feConnectomeDisplay(feSplitLoopFibers( fg ),fh,[.85 .4 .3], [], [], .35)
hold on
mctDisplayBrainSlice( niftiRead(t1File), [0 0 -10]);
set(gca,'ylim',[-126 -55],'xlim',[-5 80])
set(fh,'Position',[0 0 .45 .95],'Color',[0 0 0]);
drawnow
saveFig(fh,fullfile(saveDir,sprintf('connectome_FA_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',trackingType,lmax,bval,rep, ...
  100*diffusionModelParams(1),100*diffusionModelParams(2))));

% Display all the fiebrs for the first iteration
fh = figure('name','All fibers');
feConnectomeDisplay(feSplitLoopFibers( feGet(fe,'fibers acpc') ),fh,color, [], [], .35)
hold on
mctDisplayBrainSlice( niftiRead(t1File), [0 0 -10]);
set(gca,'ylim',[-126 -55],'xlim',[-5 80])
set(fh,'Position',[0 0 .45 .95],'Color',[0 0 0]);  drawnow
drawnow
saveFig(fh,fullfile(saveDir,sprintf('connectome_ALL_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',trackingType,lmax,bval,rep, ...
  100*diffusionModelParams(1),100*diffusionModelParams(2))));

close all
drawnow

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(h),figName));

end