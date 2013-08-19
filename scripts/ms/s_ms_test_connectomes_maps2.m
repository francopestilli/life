function s_ms_test_connectomes_maps2(trackingType,dataType,lmax,diffusionModelParams,recompute)
%
% Load FE structeres obtained by preprocessing connectomesconstrained within a
% region of interest and within the cortex and makes some basic plot of
% statistics in the connectomes
%
% See also:
%
% Franco (C) 2012 Stanford VISTA team.

if notDefined('trackingType'),trackingType = {'p'};end
if notDefined('lmax'),        lmax         = [8];end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = [1,2,3];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('cullType'),   cullType={'','culledL2'};end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_hists_maps');end

% Bins for the fiber density estimates
xBins = [1 2 4 8 16 32 64 128 256 512 1024 2048];
x     = 1:length(xBins);

% Bins for the sum of weights estimates
wxBins = [.9./(2.^[10:-1:1]) ];
wx     = 1:length(wxBins);

fontSiz    = 16;
doMaps     = 1;
doFD       = 1;
doRMSE     = 1;
figVisible = 'on';

% High-resolution Anatomy
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';
t1     = niftiRead(t1File);

for itrk = 1:length(trackingType)
    for i_lmax = 1:length(lmax)
        % Loop over Bvalues and tracking reps.
        for i_bval = 1:length(bval)
            for icull = 1:length(cullType)
                for irep = 1:length(rep)
                    % Check that the connectome were proprocessed before attempting to make a
                    % plot.
                    done = s_ms_check_processes([],trackingType{itrk},lmax(i_lmax),bval(i_bval),cullType{icull});
                    if ~all(done),
                        fprintf('\n[%s] Not all connectomes were proprocessed... not proceeding with plot.\n',mfilename);
                        return
                    end
                    
                    % File to load
                    % This is where the inputs will be loaded from
                    [feFileToLoad, fname, feLoadDir] = ...
                        msBuildFeFileName(trackingType{itrk},lmax(i_lmax),bval(i_bval),rep(irep), ...
                        diffusionModelParams,cullType{icull});
                   
                    % Rebuild the name for this special case.
                    fname =  '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax8__m_prob-500000_diffModAx100Rd0_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax8__m_prob-500000_diffModAx100Rd0_SoverS0';
                    feFileToLoad = fullfile(feLoadDir,[fname,'.mat']);

                    if (exist(feFileToLoad,'file') == 2)
                        fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)
                        load(feFileToLoad);
                    else
                        fprintf('[%s] FE file NOT found: \n%s\n ======================================== \n\n\n',mfilename,feFileToLoad)
                        keyboard
                    end
                    
                    if doFD 
                    % Get the fiber density
                    fd = feGet(fe,'fiber density');
                    end
                    
                    if doMaps
                        coords = feGet(fe,'roi coords') + 1;
                        xform  = feGet(fe,'xform img 2 acpc');
                        
                    if irep == 1
                        % Make a plot of the maps
                        slice = -70:2:-60;
                        for is = 1:length(slice)
                            if doFD
                                % Fiber density maps
                                figName = sprintf('FDM_%s_rep%i_%s_slice%i',fname,irep,cullType{icull},slice(is));
                                fh  = figure('name',figName,'visible',figVisible,'color','w');
                                img = feReplaceImageValues(nan(feGet(fe,'map size')),fd(:,icull)',coords);
                                
                                % This will be used tonormalize the fiber density plots
                                maxfd(icull,irep) = nanmax(img(:));
                                
                                % Create the nifti structure
                                ni  = niftiCreate('data',mbaNormalize(img,[0,1]) .* (maxfd(icull,irep)/maxfd(1,irep)), ...
                                    'qto_xyz',xform, ...
                                    'fname','FDM', ...
                                    'data_type',class(img));
                                sh = mbaDisplayOverlay(t1, ni, [0 slice(is) 0], [], 'hot');
                                saveMap(fh,figName,saveDir,nanmean(img(:)),nanmedian(img(:)),nanstd(img(:)),maxfd(1,irep))
                                
                                 
                                % Weigth density (sum of weights)
                                figName = sprintf('FWM_%s_rep%i_%s_slice%i',fname,irep,cullType{icull},slice(is));
                                fh  = figure('name',figName,'visible',figVisible,'color','w');
                                img = feReplaceImageValues(nan(feGet(fe,'map size')),(fd(:,3))',coords);
                                maxw(icull,irep) = nanmax(img(:));
                                minw(icull,irep) = nanmin(img(:));
                                
                                % Create the nifti structure
                                ni  = niftiCreate('data',mbaNormalize(img,[0,1]) .* (maxw(icull,irep)/maxw(1,irep)), ...
                                    'qto_xyz',xform, ...
                                    'fname','FDM', ...
                                    'data_type',class(img));
                                sh = mbaDisplayOverlay(t1, ni, [0 slice(is) 0], [], 'hot');
                                saveMap(fh,figName,saveDir,nanmean(img(:)),nanmedian(img(:)),nanstd(img(:)),100*maxw(1,irep))
                               
                            end
                            
                            if doRMSE
                                % RMSE off the model
                                figName = sprintf('RMSE_Model_%s_rep%i_%s_slice%i',fname,irep,cullType{icull},slice(is));
                                fh = figure('name',figName,'visible',figVisible,'color','w');
                                rmse = feGetRep(fe, 'vox rmse');
                                img = feReplaceImageValues(nan(feGet(fe,'map size')),rmse,coords);
                                maxr(icull,irep) = 60;
                                minr(icull,irep) = nanmin(img(:));
%                                 ni  = niftiCreate('data',mbaNormalize(img,[0,1]) .* (maxr(icull,irep)/maxr(1,irep)), ...
%                                     'qto_xyz',xform, ...
%                                     'fname','FDM', ...
%                                     'data_type',class(img));
                                ni  = niftiCreate('data',mbaNormalize(img,[0,1]), ...
                                    'qto_xyz',xform, ...
                                    'fname','FDM', ...
                                    'data_type',class(img));

                                sh = mbaDisplayOverlay(t1, ni, [0 slice(is) 0], [], 'hot');
                                saveMap(fh,figName,saveDir,nanmean(img(:)),nanmedian(img(:)),nanstd(img(:)),maxr(1,irep))
                                
                                % RMSE of the data
                                figName = sprintf('RMSE_Data_%s_rep%i_%s_slice%i',fname,irep,cullType{icull},slice(is));
                                fh = figure('name',figName,'visible',figVisible,'color','w');
                                rmse = feGetRep(fe, 'vox rmse data');
                                img = feReplaceImageValues(nan(feGet(fe,'map size')),rmse,coords);
                                maxrd(icull,irep) = 60;
                                minrd(icull,irep) = nanmin(img(:));
%                                 ni  = niftiCreate('data',mbaNormalize(img,[0,1]) .* (maxrd(icull,irep)/maxrd(1,irep)), ...
%                                     'qto_xyz',xform, ...
%                                     'fname','FDM', ...
%                                     'data_type',class(img));
                                ni  = niftiCreate('data',mbaNormalize(img,[0,1]), ...
                                    'qto_xyz',xform, ...
                                    'fname','FDM', ...
                                    'data_type',class(img));
                                sh = mbaDisplayOverlay(t1, ni, [0 slice(is) 0], [], 'hot');
                                saveMap(fh,figName,saveDir,nanmean(img(:)),nanmedian(img(:)),nanstd(img(:)),maxrd(1,irep))
                                
                                % Ratio rmse
                                figName = sprintf('Ratio_rmse_%s_rep%i_%s_slice%i',fname,irep,cullType{icull},slice(is));
                                fh = figure('name',figName,'visible',figVisible,'color','w');
                                rmse = feGetRep(fe, 'vox rmse ratio');
                                img = feReplaceImageValues(nan(feGet(fe,'map size')),rmse,coords);
                                maxrr(icull,irep) = nanmax(img(:));
                                minrr(icull,irep) = nanmin(img(:));
                                ni  = niftiCreate('data',mbaNormalize(img,[0,1]) .* (maxrr(icull,irep)/maxrr(1,irep)), ...
                                    'qto_xyz',xform, ...
                                    'fname','FDM', ...
                                    'data_type',class(img));
                                sh = mbaDisplayOverlay(t1, ni, [0 slice(is) 0], [], 'hot');
                                saveMap(fh,figName,saveDir,nanmean(img(:)),nanmedian(img(:)),nanstd(img(:)),100*maxrr(1,irep))
                            end
                            close all
                            drawnow
                        end
                    end
                    end
                    if doFD
                        % Fiber density
                        y(:,irep)     = hist(fd(:,icull),xBins);
                        ynorm(:,irep) = y(:,irep)./sum(y(:,irep));

                        % Fiber weights
                        wy(:,irep)     = hist(fd(:,3),wxBins);
                        wynorm(:,irep) = wy(:,irep)./sum(wy(:,irep));
                        
                        % Compute the dynamic range
                        dyrng{icull}(irep)  = prctile(fd(:,icull),99) / max([prctile(fd(:,icull),1),1]) ;
                        dwyrng{icull}(irep) = prctile(fd(:,3),99) / max([prctile(fd(:,3),1),0.0001]);
                    end
                end % REP REP
                
                if doFD
                    figName = sprintf('FDH_%s_%s',fname,cullType{icull});
                    fh = figure('visible',figVisible,'color','w');
                    myn{icull} = (mean(ynorm,2));
                    eyn{icull} = [myn{icull}, myn{icull}] + [-std((ynorm),[],2), std((ynorm),[],2)];
                    plot([x;x],eyn{icull}','r-','linewidth',2)
                    hold on
                    bar(x,myn{icull},'facecolor','k','edgecolor','k')
                    ylabel('Probability')
                    xlabel('Fibers per voxel')
                    
                    set(gca,'tickdir','out','ticklength',[0.025 0],'box','off','FontSize',fontSiz, ...
                        'xlim',[0 max(x)+1],'xtick',x,'xticklabel',xBins,'ylim',[0 .4])
                    saveFig(fh,fullfile(saveDir,figName),1)
                    
                    figName = sprintf('FDLL_%s_%s',fname,cullType{icull});
                    fh = figure('visible',figVisible,'color','w');
                    mllyn{icull} = log10(mean(ynorm,2));
                    ellyn{icull} = [mllyn{icull}, mllyn{icull}] + [-std(log10(ynorm),[],2), std(log10(ynorm),[],2)];
                    plot([x;x],ellyn{icull}','r-','linewidth',2)
                    hold on
                    plot(x,mllyn{icull},'ko-','markerfacecolor','k','markeredgecolor','w','markersize',8)
                    ylabel('Log-probability')
                    xlabel('Fibers per voxel')
                    set(gca,'tickdir','out','ticklength',[0.025 0],'box','off','FontSize',fontSiz,'ylim',[-3.5 0], ...
                        'xlim',[0 max(x)+1],'xtick',x,'xticklabel',xBins)
                    saveFig(fh,fullfile(saveDir,figName),1)
                    
                    
                    figName = sprintf('FWH_%s_%s',fname,cullType{icull});
                    fh = figure('visible',figVisible,'color','w');
                    wmyn{icull} = (mean(wynorm,2));
                    weyn{icull} = [wmyn{icull}, wmyn{icull}] + [-std((wynorm),[],2), std((wynorm),[],2)];
                    plot([wx;wx],weyn{icull}','r-','linewidth',2)
                    hold on
                    bar(wx,wmyn{icull},'facecolor','k','edgecolor','k')
                    ylabel('Probability')
                    xlabel('Sum of weights per voxel')
                    set(gca,'tickdir','out','ticklength',[0.025 0],'box','off','FontSize',fontSiz, ...
                        'xtick',wx,'xticklabel',wxBins,'ylim',[0 .25])
                    saveFig(fh,fullfile(saveDir,figName),1)
                    
                    figName = sprintf('FWLL_%s_%s',fname,cullType{icull});
                    fh = figure('visible',figVisible,'color','w');
                    wmllyn{icull} = log10(mean(wynorm,2));
                    wellyn{icull} = [wmllyn{icull}, wmllyn{icull}] + [-std(log10(wynorm),[],2), std(log10(wynorm),[],2)];
                    plot([wx;wx],wellyn{icull}','r-','linewidth',2)
                    hold on
                    plot(wx,wmllyn{icull},'ko-','markerfacecolor','k','markeredgecolor','w','markersize',8)
                    ylabel('Log-probability')
                    xlabel('Sum of weigths per voxel')
                    set(gca,'tickdir','out','ticklength',[0.025 0],'box','off','FontSize',fontSiz,'ylim',[-3.5 0], ...
                        'xtick',wx,'xticklabel',wxBins)
                    saveFig(fh,fullfile(saveDir,figName),1)
                end
                
            end % cull 
            if doFD
                % % FIBER DENSITY
                figName = sprintf('FDH_%s_CULLED_AND_NOT',fname);
                fh = figure('visible',figVisible,'color','w');
                colors{1} = [.6 .6 .6];
                colors{2} = [.35 .35 .35];
                for icull = 1:length(cullType)
                    plot([x;x],eyn{icull}','r-','linewidth',2)
                    hold on
                    plot(x,myn{icull},'ko-','color',colors{icull}, ...
                        'markerfacecolor',colors{icull},'markeredgecolor', ...
                        'w','markersize',8)
                    ylabel('Probability')
                    xlabel('Fibers per voxel')
                    set(gca,'tickdir','out','ticklength',[0.025 0],'box','off','FontSize',fontSiz,'ylim',[0 0.4], ...
                        'xlim',[0 max(x)+1],'xtick',x,'xticklabel',xBins);
                end
            end
            title(sprintf('Dynamic range | Before %2.2f | After %2.2f',mean(dyrng{1}),mean(dyrng{2})))
            saveFig(fh,fullfile(saveDir,figName),1)
            
            figName = sprintf('FDLL_%s_CULLED_AND_NOT',fname);
            fh = figure('visible',figVisible,'color','w');
            colors{1} = [.6 .6 .6];
            colors{2} = [.35 .35 .35];
            for icull = 1:length(cullType)
                plot([x;x],ellyn{icull}','r-','linewidth',2)
                hold on
                plot(x,mllyn{icull},'ko-','color',colors{icull}, ...
                    'markerfacecolor',colors{icull},'markeredgecolor', ...
                    'w','markersize',8)
                ylabel('Log-probability')
                xlabel('Fibers per voxel')
                set(gca,'tickdir','out','ticklength',[0.025 0],'box','off','FontSize',fontSiz,'ylim',[-3.5 0], ...
                    'xlim',[0 max(x)+1],'xtick',x,'xticklabel',xBins)
            end
            saveFig(fh,fullfile(saveDir,figName),1)

            % % FIBER WEIGHTS
            figName = sprintf('FWH_%s_CULLED_AND_NOT',fname);
            fh = figure('visible',figVisible,'color','w');
            colors{1} = [.6 .6 .6];
            colors{2} = [.35 .35 .35];
            for icull = 1:length(cullType)
                plot([wx;wx],weyn{icull}','r-','linewidth',2)
                hold on
                plot(wx,wmyn{icull},'ko-','color',colors{icull}, ...
                    'markerfacecolor',colors{icull},'markeredgecolor', ...
                    'w','markersize',8)
                ylabel('Probability')
                xlabel('Sum of weights in voxel')
                set(gca,'tickdir','out','ticklength',[0.025 0],'box','off','FontSize',fontSiz,'ylim',[0 0.4], ...
                    'xtick',wx,'xticklabel',wxBins,'ylim',[0 .25])
            % Compute the dynamic range:
            maxx = wmyn{icull}>0;
            dyrng{icull} = diff(minmax(wx(maxx)));
            end
            title(sprintf('Dynamic range | Before %2.2f | After %2.2f',mean(dwyrng{1}),mean(dwyrng{2})))
     
            saveFig(fh,fullfile(saveDir,figName),1)
            
            figName = sprintf('FWLL_%s_CULLED_AND_NOT',fname);
            fh = figure('visible',figVisible,'color','w');
            colors{1} = [.6 .6 .6];
            colors{2} = [.35 .35 .35];
            for icull = 1:length(cullType)
                plot([wx;wx],wellyn{icull}','r-','linewidth',2)
                hold on
                plot(wx,wmllyn{icull},'ko-','color',colors{icull}, ...
                    'markerfacecolor',colors{icull},'markeredgecolor', ...
                    'w','markersize',8)
                ylabel('Log-probability')
                xlabel('Sum of weights in voxel')
                set(gca,'tickdir','out','ticklength',[0.025 0],'box','off','FontSize',fontSiz,'ylim',[-3.5 0], ...
                    'xtick',wx,'xticklabel',wxBins)
                
            end
            saveFig(fh,fullfile(saveDir,figName),1)
            %close all
            drawnow
            end
        end     % bval
    end         % lmax
end             % tracking type
end             % Main function


%---------------------------------%
function saveMap(fh,figName,saveDir,M,m,SD,maxfd)
% This helper function saves two figures for each map and eps with onlythe
% axis and a jpg with only the brain slice.
% The two can then be combined in illustrator.
%
% First we save only the slice as jpeg.
set(gca,'fontsize',16,'ztick',[-20 -10 0 10 20], ...
    'xtick',[0 10 20 30 40 50], ...
    'xlim',[-5 55],'zlim',[-22 22],'tickdir','out','ticklength',[0.025 0])
axis off
saveFig(fh,fullfile(saveDir,'maps',figName),'tiff')
saveFig(fh,fullfile(saveDir,'maps',figName),'png')

% Then we save the slice with the axis as
% eps. This will only generate the axis
% that can be then combined in illustrator.
axis on
grid off

title(sprintf('mean %2.2f | median %2.2f | SD %2.2f', ...
    M,m,SD),'fontsize',16)
zlabel('Z (mm)','fontsize',16)
xlabel('X (mm)','fontsize',16)
cmap = colormap(hot(255));
colorbar('ytick',linspace(0,1,5),'yticklabel', ...
    {1, num2str(ceil(maxfd/8)), num2str(ceil(maxfd/4)), ...
    num2str(ceil(maxfd/2)), num2str(ceil(maxfd))}, ...
    'tickdir','out','ticklength',[0.025 0],'fontsize',16)
saveFig(fh,fullfile(saveDir,'maps',figName),1)
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