function s_ms_plot_culling_results_across_iterations(fitType,baseDir)
%
% This function makes some plots of the results of a culling of connectomes
% across iterations.
%
% It is in dev state, it might improve if these plots become useful.
%
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

close all
clear h

if notDefined('fitType'), fitType = 'L2';end
os  = {sprintf('oTensor%s',fitType),sprintf('oDet%s',fitType),sprintf('oProb%s',fitType)};
fes = {sprintf('feTensor%s',fitType),sprintf('feDet%s',fitType),sprintf('feProb%s',fitType)};

% Load the results from disk:
for ii = 1:length(fes)
    if ~exist(fes{ii},'var')
        disp('loading culling results...')
        load(fullfile('/home/frk/',[fes{ii}(1:end-2),'.mat']))
    end
end

% Data directories:
% used copy: /home/frk
if notDefined('baseDir')
    baseDir = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/culling_results_l1_l2/';
end

% Set up plots and conditions
fontSiz = 18;
%colors = {[.125 .5 .5],[.5 .5 .125],[.5 .125 .5]};
colors = {[.5 .4 .65],[.6 .45 .4],[.45 .6 .35]};
r2 = cell(length(os),1);rmse=r2;rrmse=r2;
nFibs = r2;

% Compute several statistics, plot some of them.
for ii=1:length(os)
    useIdx =  eval(sprintf('~isnan(%s.r2);',os{ii}));
    r2{ii} = 100*eval(sprintf('%s.r2xv( useIdx );',os{ii}));    
    rmse{ii} = eval(sprintf('%s.rmsexv( useIdx );',os{ii}));
    rrmse{ii} = eval(sprintf('%s.rrmse( useIdx );',os{ii}));
    nFibs{ii} = eval(sprintf('%s.numFibers( useIdx );',os{ii}));
    drmse{ii} = median(feGetRep(eval(sprintf('%s',fes{ii})),'vox rmse data'));
    %mrmse{ii} = median(feGetRep(eval(sprintf('%s',fes{ii})),'vox rmse'));
    %mr2{ii}   = median(feGetRep(eval(sprintf('%s',fes{ii})),'vox r2'));
    dsd{ii}     = var(feGetRep(eval(sprintf('%s',fes{ii})),'dsigdemeaned'));
    dsig2{ii}  = feGetRep(eval(sprintf('%s',fes{ii})),'dsigdemeaned');
    dsig1{ii}  = feGet(eval(sprintf('%s',fes{ii})),'dsigdemeaned');
    dr2{ii}   = median(feGetRep(eval(sprintf('%s',fes{ii})),'vox r2 data'));
    mmpVx = feGet(eval(sprintf('%s',fes{ii})),'xform img 2 acpc');
    mmpVx = mmpVx(1:3,1:3);
    mmpVx = diag(mmpVx);
    vxNum{ii} = size(unique(floor(mrAnatXformCoords(feGet(eval(sprintf('%s',fes{ii})), ...
        'xform acpc 2 img'), ...
        feGet(eval(sprintf('%s',fes{ii})),'roi coords'))),'rows'),1);
    totWMvol{ii} = vxNum{ii}*prod(mmpVx);% in mm^3
end

for ii=1:length(os)
    fname = sprintf('Data_dSig_reliability_%s_%s',fitType,fes{ii});
    fh = mrvNewGraphWin(fname);
    plot(dsig1{ii},dsig2{ii},'o','MarkerFaceColor',colors{ii},'MarkerEdgeColor',colors{ii});
    axis square
    set(gca,'ylim',[-300 300],'ytick',[-300 -150 0 150 300], ...
        'ylim',[-300 300],'ytick',[-300 -150 0 150 300], ...
        'tickdir','out','box','off','FontSize',fontSiz)
    ylabel('Demeaned diffusion signal data 1');
    xlabel('Demeaned diffusion signal data 2');
    saveFig(fh,fullfile(baseDir,fname)) 
end

fname = sprintf('Data_dSig_hist_%s_%s',fitType,fes{ii});
fh = mrvNewGraphWin(fname);
for ii=length(os):-1:1
    [y,x] = hist(dsig1{ii},120);
    hh = bar(x,y./totWMvol{ii},'facecolor',colors{ii},'edgecolor','w','linewidth',.01);
    set(get(hh,'Children'),'FaceAlpha',1,'EdgeAlpha',1)
    hold on
end
set(gca,  'tickdir','out','box','off','FontSize',fontSiz)
ylabel('Percent white-matter in connectome');
xlabel('Demeaned diffusion signal');
saveFig(fh,fullfile(baseDir,fname))


fname = sprintf('Data_R2_data_reliability_%s',fitType);
fh = mrvNewGraphWin(fname);
for ii=1:length(os)
    bar(ii,100*dr2{ii},'FaceColor',colors{ii});
    hold on
end
set(gca,'ylim',[20 50],'ytick',[20 30 40 50], ...
    'xtick',[1 2 3],'xticklabel',[os],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('Percent variance explained from data 1 to data 2');xlabel('Algorithm type')
saveFig(fh,fullfile(baseDir,fname)) 

fname = ('Variance_in_the_data');
fh = mrvNewGraphWin(fname);
for ii=1:length(os)
    bar(ii,dsd{ii},'FaceColor',colors{ii});
    hold on
end
set(gca,'ylim',[2000 2500],'ytick',[2000 2250 2500], ...
    'xtick',[1 2 3],'xticklabel',[os],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('Demeaned diffusion signal variance');xlabel('Algorithm type')
saveFig(fh,fullfile(baseDir,fname))

fname = sprintf('Data_RMSE_data_reliability)_%s',fitType);
fh = mrvNewGraphWin(fname);
for ii=1:length(os)
    bar(ii,drmse{ii},'FaceColor',colors{ii});
    hold on
end
set(gca,'ylim',[28 30],'ytick',[28 29 30], ...
    'xtick',[1 2 3],'xticklabel',[os],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('RMSE');xlabel('Algorithm type')
saveFig(fh,fullfile(baseDir,fname)) 

fname = sprintf('R2_%s',fitType);
fh = mrvNewGraphWin(fname);
for ii=1:length(os)
    h(ii) = semilogx(1:length(r2{ii}),r2{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii})    ;
    hold on

end
set(gca,'ylim',[20 35],'ytick',[20 25 30 35],'xtick',[0 10 100 1000],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('Percent variance explained');xlabel('iteration number')
saveFig(fh,fullfile(baseDir,fname)) 

fname = sprintf('rmse_%s',fitType);
fh = mrvNewGraphWin(fname);
for ii=1:length(os)
    h(ii) = semilogx(1:length(rmse{ii}),rmse{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii})   ;
    hold on
end
set(gca,'ylim',[26 32],'ytick',[26 30 32],'xtick',[0 10 100 1000],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('RMSE');xlabel('iteration number')
legend(h,os,'box','off')
saveFig(fh,fullfile(baseDir,fname)) 

fname = sprintf('Rrmse_%s',fitType);
fh = mrvNewGraphWin(fname);
for ii=1:length(os)
    h(ii) = semilogx(1:length(rrmse{ii}),rrmse{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii})   ;
    hold on
end
set(gca,'ylim',[0.75 1.25],'ytick',[.75 1 1.25],'xtick',[0 10 100 1000],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('R_{rmse}');xlabel('iteration number')
legend(h,os,'box','off')
saveFig(fh,fullfile(baseDir,fname)) 

fname = sprintf('Deleted_fibers_%s',fitType);
fh = mrvNewGraphWin(fname);
for ii=1:length(os)
    h(ii) = loglog(.0001+(1:length(nFibs{ii})),nFibs{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii});
    hold on
end
set(gca,'xlim',[1 1000],'ylim',[1000 64000],'ytick',[1000 2000 4000 8000 16000 32000 64000], ...
    'xtick',[1 10 100 1000],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('Number of fascicles in connectome');xlabel('iteration number')
legend(h,os,'box','off')
saveFig(fh,fullfile(baseDir,fname)) 

fname = sprintf('Deleted_fibers_versus_rmse_%s',fitType);
fh = mrvNewGraphWin(fname);
for ii=1:length(os)
    h(ii) = semilogx(nFibs{ii},rmse{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii}) ;
    hold on
end
set(gca,'ylim',[26 32],'ytick',[26 29 32],'xlim',[1000 64000],'xtick',[1000 2000 4000 8000 16000 32000 64000],'tickdir','out','box','off','FontSize',fontSiz)
xlabel('Number of fascicles in connectome');ylabel('rmse')
legend(h,os,'box','off')
saveFig(fh,fullfile(baseDir,fname)) 

fname = sprintf('Deleted_fibers_versus_R2_%s',fitType);
fh = mrvNewGraphWin(fname);
for ii=1:length(os)
    h(ii) = semilogx(nFibs{ii},r2{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii})   ;
    hold on
end
set(gca,'ylim',[20 40],'ytick',[20 30 40],'xlim',[1000 64000],'xtick',[1000 2000 4000 8000 16000 32000 64000],'tickdir','out','box','off','FontSize',fontSiz)
xlabel('Number of fascicles in connectome');ylabel('R^2')
legend(h,os,'box','off')
saveFig(fh,fullfile(baseDir,fname)) 

fname = sprintf('Rrmse_vs_R2_%s',fitType);
fh = mrvNewGraphWin(fname);
for ii=1:length(os)
    h(ii) = semilogx(rrmse{ii},r2{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii})   ;
    hold on
end
set(gca,'ylim',[20 40],'ytick',[20 30 40],'xlim',[0.9 1.1],'xtick',[0.95 1 1.05],'tickdir','out','box','off','FontSize',fontSiz)
xlabel('R_{rmse}');ylabel('R^2')
legend(h,os,'box','off')
saveFig(fh,fullfile(baseDir,fname)) 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)

printCommand = sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName);
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

% do the printing here:
eval(printCommand);

end