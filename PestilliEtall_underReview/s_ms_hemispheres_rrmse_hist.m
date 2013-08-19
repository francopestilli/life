function s_ms_hemispheres_rrmse_hist(trackingType,lmax,bval,rep,volume)
%
% s_ms_hemispheres_rrmse_hist(trackingType,lmax,bval,rep)
%
% Load a series of FE structures for conectomes of the occipital lobe.
% Saves figures of the occipital diffusion signal (measured, predicted, error).
%
% This is part of a series of reproducible science scritps to be published
% with the LiFE mansucript.
%
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

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
dsig_colormap = 'bone';
figVisible    = 'on';
% Make a histogram plot of the fiber density and of the Rrmse
nBins= logspace(log10(.5),log10(2),25);
x     = 1:length(nBins);
% Directory to save the figures
saveDirF = fullfile(saveDir,'rrmse_hist');

% High-resolution Anatomy
t1File    = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';
t1        = niftiRead(t1File);
basedir   = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/';
feType{1} = 'fe_hemispheres/fe_culled_FP_150_B2000_LMAX8_right.mat';
feType{2} = 'fe_hemispheres/fe_culled_FP_150_B2000_LMAX8_left.mat';

for irep = 2:3
    rrmse = [];
    % Information on the path to the files to load.
    feFileToLoad{1} = fullfile(basedir,sprintf('life_mrtrix_rep%i',irep),feType{1});
    feFileToLoad{2} = fullfile(basedir,sprintf('life_mrtrix_rep%i',irep),feType{2});
    
    for ih = 1:2
        disp('Loading the FE structure...')
        load(feFileToLoad{ih});
        
        % Rrmse
        rSig{ih}    = feGetRep(fe,'voxrmseratio');
    end
    
    rrmse =  rSig{:};
    clear rSig
    
    % Build a histrogram
    [y(:,irep-1),x] = hist(rrmse,nBins);
    nSum = sum(y(:,irep-1));
    y(:,irep-1) = y(:,irep-1)./nSum;
end


% Compute mean and std
my = mean(y,2);
yerr = [my,my] + 3*[-std(y,[],2),std(y,[],2)];
my = my';
yerr = yerr';

figName = sprintf('Rrmse_hist');
fh = figure('name',figName,'visible',figVisible,'color','w',...
            'Units','normalized','Position',[0 .1  0.35  0.95]);

colors{1} = [.35 .35 .35];
px =  x( x <= 1 );
px = [px px( end )];
py = [my( x <= 1 ) 0];
pp = patch(px,py,[.8 .8 .8],'edgecolor',[.8 .8 .8]);
hold on
plot([1 1],[0 .16],'k--','linewidth',2)
plot(x,my','o-','color',colors{1}, ...
    'markerfacecolor', colors{1}, ...
    'markeredgecolor','w',...
    'markersize',18,'linewidth',2)
plot([x;x],yerr,'r-','linewidth',3)
set(gca,'tickdir','out','box','off', ...
    'fontsize',20,'ylim',[0 .16],'ytick',[0 .08 .16], ...
    'xtick',[.5 1 2],  'xscale','log','ticklength',[0.025 0])
ylabel('Probability','fontsize',20)
xlabel('R_{rmse}','fontsize',20)
title(sprintf('Proportion R_{rmse}<= 1: %2.3f',sum(y( x <= 1 ))),'fontsize',20)
saveFig(fh,fullfile(saveDirF,figName),1)

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