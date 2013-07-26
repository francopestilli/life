function s_ms_test_connectomes_display_dwi(trackingType,lmax,bval,rep,volume)
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
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_hemispheres_signal');end
doFD = 0;

% Figures, slices, axis and plotting.
dirs   = [65];
slices = [-70];
xlim   = [-70 67];
zlim   = [-22 75];
dsig_colormap = 'winter';
figVisible    = 'on';

% High-resolution Anatomy
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';
t1     = niftiRead(t1File);

% Information on the path to the files to load.
feFileToLoad{1} = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/fe_hemispheres/fe_culled_FP_150_B2000_LMAX8_right.mat';
feFileToLoad{2} = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/fe_hemispheres/fe_culled_FP_150_B2000_LMAX8_left.mat';
for isl = 1:length(slices)
    for idir = 1:length(dirs)
        volume = [0 slices(isl) 0 dirs(idir)];
        for ih = 1:2
            disp('Loading the FE structure...')
            load(feFileToLoad{ih});
            
            % Indices into the 4D volume toaddress the signal without the
            % B0 measurements
            volSiz = feGet(fe,'volumesize')-[0 0 0 10];

            % Get the xform and the coordinates
            xform   = feGet(fe,'xform img 2 acpc');
            coords{ih} = feGet(fe,'roicoords')+1; % This is weird, it appears that i need to add 1 to all the coordinates
             
            if doFD
                % Get the fiber density
                fd{ih} = feGet(fe,'fiber density');
                fdImg{ih} = feValues2volume(fd{ih}(:,1)',coords{ih},volSiz);
            end
            
            % Get the signal into an image.
            dSig{ih}     = feGet(fe,'dsigdemeanedvox');
            dSigImg1{ih} = feValues2volume(dSig{ih},coords{ih},volSiz);
            
            % Get the signal into an image.
            dSig{ih}     = feGetRep(fe,'dsigdemeanedvox');
            dSigImg2{ih} = feValues2volume(dSig{ih},coords{ih},volSiz);
            
            % Predicted signal
            pSig{ih}    = feGet(fe,'psigfvox');
            pSigImg{ih} = feValues2volume(pSig{ih},coords{ih},volSiz);
            
            % The last terms do not have a 4th dimension           
            volSiz  = volSiz(1:3);

            % Error model
            eSig{ih}    = feGetRep(fe,'voxrmse');
            eSigImg{ih} = feValues2volume(eSig{ih},coords{ih},volSiz);
            
            % Error data
            edSig{ih}    = feGetRep(fe,'voxrmsedata');
            edSigImg{ih} = feValues2volume(edSig{ih},coords{ih},volSiz);
            
            % Rrmse
            rSig{ih}    = feGetRep(fe,'voxrmseratio');
            rSigImg{ih} = feValues2volume(rSig{ih},coords{ih},volSiz);
        end
        if doFD
            % Combine the informantion of the left and right hemisphere
            fdImgB = nan(size(fdImg{1}));
            fdImgB(~isnan(fdImg{1})) = fdImg{1}(~isnan(fdImg{1}));
            fdImgB(~isnan(fdImg{2})) = fdImg{2}(~isnan(fdImg{2}));
        end
        
        % Combine the informantion of the left and right hemisphere
        dSigImgB1 = nan(size(dSigImg1{1}));
        dSigImgB1(~isnan(dSigImg1{1})) = dSigImg1{1}(~isnan(dSigImg1{1}));
        dSigImgB1(~isnan(dSigImg1{2})) = dSigImg1{2}(~isnan(dSigImg1{2}));
        
        % Combine the informantion of the left and right hemisphere
        dSigImgB2 = nan(size(dSigImg2{1}));
        dSigImgB2(~isnan(dSigImg2{1})) = dSigImg2{1}(~isnan(dSigImg2{1}));
        dSigImgB2(~isnan(dSigImg2{2})) = dSigImg2{2}(~isnan(dSigImg2{2}));
        
        % Combine the informantion of the left and right hemisphere
        pSigImgB = nan(size(pSigImg{1}));
        pSigImgB(~isnan(pSigImg{1})) = pSigImg{1}(~isnan(pSigImg{1}));
        pSigImgB(~isnan(pSigImg{2})) = pSigImg{2}(~isnan(pSigImg{2}));
        
        % Combine the informantion of the left and right hemisphere
        eSigImgB = nan(size(eSigImg{1}));
        eSigImgB(~isnan(eSigImg{1})) = eSigImg{1}(~isnan(eSigImg{1}));
        eSigImgB(~isnan(eSigImg{2})) = eSigImg{2}(~isnan(eSigImg{2}));
         
        % Combine the informantion of the left and right hemisphere
        edSigImgB = nan(size(edSigImg{1}));
        edSigImgB(~isnan(edSigImg{1})) = edSigImg{1}(~isnan(edSigImg{1}));
        edSigImgB(~isnan(edSigImg{2})) = edSigImg{2}(~isnan(edSigImg{2}));
        
        % Combine the informantion of the left and right hemisphere
        rSigImgB = nan(size(rSigImg{1}));
        rSigImgB(~isnan(rSigImg{1})) = rSigImg{1}(~isnan(rSigImg{1}));
        rSigImgB(~isnan(rSigImg{2})) = rSigImg{2}(~isnan(rSigImg{2}));
        
        if doFD
            % Create the nifti structure
            niFD  = niftiCreate('data',fdImgB(:,:,:,volume(end)), ...
                'qto_xyz',xform, ...
                'fname','fd', ...
                'data_type',class(fdImgB));
        end
        
        % Create the nifti structure
        niMs1  = niftiCreate('data',dSigImgB1(:,:,:,volume(end)), ...
            'qto_xyz',xform, ...
            'fname','dwi_measured_signal', ...
            'data_type',class(dSigImgB1));
        
        % Create the nifti structure
        niMs2  = niftiCreate('data',dSigImgB2(:,:,:,volume(end)), ...
            'qto_xyz',xform, ...
            'fname','dwi_measured_signal', ...
            'data_type',class(dSigImgB2));
        
        % Create the nifti structure
        niP  = niftiCreate('data',pSigImgB(:,:,:,volume(end)), ...
            'qto_xyz',xform, ...
            'fname','dwi_predicted_signal', ...
            'data_type',class(pSigImgB));
        
        % Create the nifti structure
        niE  = niftiCreate('data',eSigImgB, ...
            'qto_xyz',xform, ...
            'fname','dwi_rmse_signal', ...
            'data_type',class(eSigImgB));
        
        % Create the nifti structure
        niED  = niftiCreate('data',edSigImgB, ...
            'qto_xyz',xform, ...
            'fname','dwi_rmse_signal', ...
            'data_type',class(edSigImgB));
        
        % Create the nifti structure
        niR  = niftiCreate('data',rSigImgB, ...
            'qto_xyz',xform, ...
            'fname','dwi_rmse_signal', ...
            'data_type',class(rSigImgB));
        
        % Directory to save the figures
        saveDirF = fullfile(saveDir,['slice',num2str(slices(isl))],['dir',num2str(dirs(idir))]);
        
        if doFD
        % Measured signal 1
        figName = sprintf('FiberDensity_slice%i_dir%i_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',...
                  slices(isl),dirs(idir),trackingType,lmax,bval,rep, ...
                  100*diffusionModelParams(1),100*diffusionModelParams(2));
        fh = figure('name',figName,'visible',figVisible,'color','w');
        sh = mbaDisplayOverlay(t1, niFD, volume(1:3), [], 'hot');
                    
        % Tick marks for the colorbar
        mm = minmax(niFD.data(:));
        barticks = [mm(1) 0 mm(2)];
                       
        % Information to display in the title
        M  = nanmean(  niFD.data(:));
        m  = nanmedian(niFD.data(:));
        SD = nanstd(   niFD.data(:));
        saveMap(fh, figName, saveDirF,M,m,SD,barticks,xlim,zlim,'hot');
        end
        
        % Measured signal 1
        figName = sprintf('dwi_measured_signal_1_slice%i_dir%i_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',...
                  slices(isl),dirs(idir),trackingType,lmax,bval,rep, ...
                  100*diffusionModelParams(1),100*diffusionModelParams(2));
        fh = figure('name',figName,'visible',figVisible,'color','w');
        sh = mbaDisplayOverlay(t1, niMs1, volume(1:3), [], dsig_colormap);
                    
        % Tick marks for the colorbar
        mm = minmax(niMs1.data(:));
        barticks = [mm(1) 0 mm(2)];
                       
        % Information to display in the title
        M  = nanmean(  niMs1.data(:));
        m  = nanmedian(niMs1.data(:));
        SD = nanstd(   niMs1.data(:));
        saveMap(fh, figName, saveDirF,M,m,SD,barticks,xlim,zlim,dsig_colormap);
        
        % Measured signal 2
        figName = sprintf('dwi_measured_signal_2_slice%i_dir%i_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',...
                  slices(isl),dirs(idir),trackingType,lmax,bval,rep, ...
                  100*diffusionModelParams(1),100*diffusionModelParams(2));
        fh = figure('name',figName,'visible',figVisible,'color','w');
        sh = mbaDisplayOverlay(t1, niMs2, volume(1:3), [], dsig_colormap);
          
        % Tick marks for the colorbar
        mm = minmax(niMs2.data(:));
        barticks = [mm(1) 0 mm(2)];%[-200 0 200];
                      
        % Information to display in the title
        M  = nanmean(  niMs2.data(:));
        m  = nanmedian(niMs2.data(:));
        SD = nanstd(   niMs2.data(:));
        saveMap(fh, figName, saveDirF,M,m,SD,barticks,xlim,zlim,dsig_colormap);
        
        % Predicted signal
        figName = sprintf('dwi_predicted_signal_slice%i_dir%i_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',...
                  slices(isl),dirs(idir),trackingType,lmax,bval,rep, ...
                  100*diffusionModelParams(1),100*diffusionModelParams(2));
        fh = figure('name',figName,'visible',figVisible,'color','w');
        sh = mbaDisplayOverlay(t1, niP, volume(1:3),[],dsig_colormap);
        
        % Tick marks for the colorbar
        mm = minmax(niP.data(:));
        barticks = [mm(1) 0 mm(2)];%[-200 0 200];
        
        % Information to display in the title
        M  = nanmean(  niP.data(:));
        m  = nanmedian(niP.data(:));
        SD = nanstd(   niP.data(:));
        saveMap(fh, figName, saveDirF,M,m,SD,barticks,xlim,zlim,dsig_colormap);
    end
         
    if doFD
        saveDirFD = fullfile(saveDir,'fiber_density');
        % Make a histogram plot of the fiber density and of the Rrmse
        xBins = [0 1 2 4 8 16 32 64];
        x     = 1:length(xBins);
        figName = sprintf('FiberDensity_hist');
        fh = figure('name',figName,'visible',figVisible,'color','w');
        colors{1} = [.6 .6 .6];
        colors{2} = [.35 .35 .35];
        
        % Fiber density
        FD = fd{:};
        y     = hist(FD(:,1),xBins);
        ynorm = y./sum(y);
        
        % Compute the dynamic range
        dyrng  = prctile(FD(:,1),99) / max([prctile(FD(:,1),1),1]) ;
        plot(x,ynorm,'ko-','color',colors{2}, ...
            'markerfacecolor',colors{2},'markeredgecolor', ...
            'w','markersize',8)
        ylabel('Probability')
        xlabel('Fibers per voxel')
        set(gca,'tickdir','out','ticklength',[0.025 0],'box','off','FontSize',20,'ylim',[0 0.4], ...
            'xlim',[0 max(x)+1],'xtick',x,'xticklabel',xBins);
        title(sprintf('Dynamic range %2.2f',dyrng))
        saveFig(fh,fullfile(saveDirFD,figName),1)
    end
    
    % Directory to save the figures
    saveDirF = fullfile(saveDir,'errors');
    
    % Make a histogram plot of the fiber density and of the Rrmse
    nBins= logspace(log10(.5),log10(2),25);   
    x     = 1:length(nBins);
    Rrmse = rSig{:};
    [y,x] = hist(Rrmse,nBins);
    nSum = sum(y);
    y = y./nSum;
    
    figName = sprintf('Rrmse_hist');
    fh = figure('name',figName,'visible',figVisible,'color','w');
    colors{1} = [.35 .35 .35];
    px = x(x<=1);
    px =[px px(end)];
    py = [y(x<=1) 0];
    pp = patch(px,py,[.8 .8 .8],'edgecolor',[.8 .8 .8]);
    hold on
    plot([1 1],[0 .16],'k--')
    plot(x,y,'o-','color',colors{1}, ...
        'markerfacecolor', colors{1}, ...
        'markeredgecolor','w',...
        'markersize',18,'linewidth',2)
    set(gca,'tickdir','out','box','off', ...
        'fontsize',20,'ylim',[0 .16],'ytick',[0 .08 .16], ...
        'xtick',[.5 1 2],  'xscale','log','ticklength',[0.025 0])
    ylabel('Probability','fontsize',20)
    xlabel('R_{rmse}','fontsize',20)
    title(sprintf('Proportion R_{rmse}<= 1: %2.3f',sum(y(x<=1))),'fontsize',20)
    saveFig(fh,fullfile(saveDirF,figName),1)
        
        % RMSE
        figName = sprintf('dwi_rmse_model_slice%i_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',....
            slices(isl),trackingType,lmax,bval,rep, ...
            100*diffusionModelParams(1),100*diffusionModelParams(2));
        fh = figure('name',figName,'visible',figVisible,'color','w');
        sh = mbaDisplayOverlay(t1, niE, volume(1:3));
        
        % Tick marks for the colorbar
        mm = round(minmax(niE.data(:)));
        barticks = [mm(1) mean(mm) mm(2)];
        
        % Information to display in the title
        M  = nanmean(  niE.data(:));
        m  = nanmedian(niE.data(:));
        SD = nanstd(   niE.data(:));
        saveMap(fh, figName, saveDirF,M,m,SD,barticks,xlim,zlim,'hot');
         
        % RMSE data
        figName = sprintf('dwi_rmse_data_slice%i_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',....
            slices(isl),trackingType,lmax,bval,rep, ...
            100*diffusionModelParams(1),100*diffusionModelParams(2));
        fh = figure('name',figName,'visible',figVisible,'color','w');
        sh = mbaDisplayOverlay(t1, niED, volume(1:3));
        
        % Tick marks for the colorbar
        mm = round(minmax(niED.data(:)));
        barticks = [mm(1) mean(mm) mm(2)];
        
        % Information to display in the title
        M  = nanmean(  niED.data(:));
        m  = nanmedian(niED.data(:));
        SD = nanstd(   niED.data(:));
        saveMap(fh, figName, saveDirF,M,m,SD,barticks,xlim,zlim,'hot');
        
        % Rrmse
        figName = sprintf('dwi_Rrmse_slice%i_%s_lmax%i_bval%i_rep%i_diffMode%i_%i',...
                  slices(isl),trackingType,lmax,bval,rep, ...
                  100*diffusionModelParams(1),100*diffusionModelParams(2));
        fh = figure('name',figName,'visible',figVisible,'color','w');
        sh = mbaDisplayOverlay(t1, niR, volume(1:3));
        
        % Tick marks for the colorbar
        mm = minmax(niR.data(:));
        barticks = [0.5 1 2];
        
        % Information to display in the title
        M  = nanmean(  niR.data(:));
        m  = nanmedian(niR.data(:));
        SD = nanstd(   niR.data(:));
        saveMap(fh, figName, saveDirF,M,m,SD,barticks,xlim,zlim,'hot');
    
 drawnow
end

end

%---------------------------------%
function saveMap(fh,figName,saveDir,M,m,SD,barticks,xlim,zlim,mapType)
% This helper function saves two figures for each map and eps with onlythe
% axis and a jpg with only the brain slice.
% The two can then be combined in illustrator.
%
% We save only the slice as jpeg.
set(fh,'Units','normalized','Position',[0 .1  0.35  0.95]);
set(gca,'fontsize',20,'ztick',[-40 -20 0 20 40 60], ...
    'xtick',[-50 -25 0 25 50], ...
    'xlim',xlim,'zlim',zlim,'tickdir','out','ticklength',[0.025 0])
axis off
saveFig(fh,fullfile(saveDir,figName),'png')

% Then we save the slice with the axis as
% eps. This will only generate the axis
% that can be then combined in illustrator.
axis on
grid off

% Title and lables information
title(sprintf('mean %2.2f | median %2.2f | SD %2.2f', M,m,SD),'fontsize',16)
zlabel('Z (mm)','fontsize',20)
xlabel('X (mm)','fontsize',20)

% Build a colormap
cmap = colormap(eval(sprintf('%s(255)',mapType)));
ch = colorbar('ytick',linspace(0,1,3),'yticklabel', ...
    barticks, 'tickdir','out', ...
    'ticklength', [0.025 0], 'fontsize',20);
drawnow
saveFig(fh,fullfile(saveDir,figName),'eps');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName,eps)
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

switch eps
    case {0,'jpeg'}
        eval(sprintf('print(%s, ''-djpeg90'', ''-opengl'', ''%s'')', num2str(h),figName));
    case {1,'eps'}
        eval(sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName));
    case 'png'
        eval(sprintf('print(%s, ''-dpng'',''-r500'', ''%s'')', num2str(h),figName));
    case 'tiff'
        eval(sprintf('print(%s, ''-dtiff'',''-r500'', ''%s'')', num2str(h),figName));
    case 'bmp'
        eval(sprintf('print(%s, ''-dbmp256'',''-r500'', ''%s'')', num2str(h),figName));
  otherwise
end

end