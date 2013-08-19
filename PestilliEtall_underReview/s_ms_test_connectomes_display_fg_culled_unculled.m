function s_ms_test_connectomes_display_fg_culled_unculled(trackingType,lmax,bval,rep)
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


if notDefined('trackingType'),trackingType = 'p';end
if notDefined('lmax'),        lmax         = 10;end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = 1;end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_culling_figs');end

% Perecentile cut-off for the weights
color = [.8 .7 .88];

% High-resolution Anatomy
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';

culledType = {'culledL2',''};

for iCul = 1:length(culledType)
% Information on the path to the files to load.
% This is where the inputs will be loaded from
[feFileToLoad, fname] = msBuildFeFileName(trackingType,lmax,bval,rep,diffusionModelParams,culledType{iCul});

% Get the fe structure
fprintf('[%s] Loading %s...\n',mfilename, feFileToLoad);
load(feFileToLoad);

% We will now show and save out three connectomes. The full one, the one
% with the best fibers and the one with only the fibers rejected by the
% life fit.

% Find the fibers that have non-zero weight for life (Hits)
% Rough way to the the voxel-wise weights. Think more what to do wit this.
w.all = (fe.life.fit.weights);

% Find the positive non-zero weights
w.good  = w.all > 0 ;

% Extract the good fibers
fg = fgExtract(feGet(fe,'fibers acpc'),find(w.good),'keep');

% Dispaly the good fibers 
fh = figure('name','Accepted fibers');
feConnectomeDisplay(feSplitLoopFibers( fg ),fh,color, [], [], .235);
hold on
mctDisplayBrainSlice( niftiRead(t1File), [0 0 -10]);
set(gca,'ylim',[-110 -59.5],'xlim',[-2 52])
set(fh,'Position',[0 0 .45 .95],'Color',[0 0 0]);
drawnow
saveFig(fh,fullfile(saveDir,sprintf('connectome_%s_%s',fname,culledType{iCul})));

set(gca,'ylim',[-97.8 -80],'xlim',[.75 38])
set(fh,'Position',[0 0 .45 .95],'Color',[0 0 0]);
drawnow
saveFig(fh,fullfile(saveDir,sprintf('connectome_%s_%s_ZOOM',fname,culledType{iCul})));
close all
drawnow
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end

eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(h),figName));

end