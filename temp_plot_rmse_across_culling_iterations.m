close all
clear h
% colors = {[.125 .5 .5],[.125 .5 .3],[.5 .5 .125],[.5 .3 .125],[.5 .125 .5],[.3 .125 .5]};
% os = {'oTensorL1','oTensorL2','oDetL1','oDetL2','oProbL1','oProbL2'};

% Data directories:
% backup: /azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/culling_results_l1_l2/
% used copy: /home/frk

fontSiz = 18;
colors = {[.125 .5 .5],[.5 .5 .125],[.5 .125 .5]};
os = {'oTensorL1','oDetL1','oProbL1'};
fes = {'feTensorL1','feDetL1','feProbL1'};
r2 = cell(length(os),1);rmse=r2;rrmse=r2;
nFibs = r2;

% Make a plot of the the cross-validate R2
for ii=1:length(os)
    useIdx =  eval(sprintf('~isnan(%s.r2);',os{ii}));
    
    r2{ii} = 100*eval(sprintf('%s.r2xv( useIdx );',os{ii}));    
    rmse{ii} = eval(sprintf('%s.rmsexv( useIdx );',os{ii}));
    rrmse{ii} = eval(sprintf('%s.rrmse( useIdx );',os{ii}));
    nFibs{ii} = eval(sprintf('%s.numFibers( useIdx );',os{ii}));
    drmse{ii} = median(feGetRep(eval(sprintf('%s',fes{ii})),'vox rmse data'));
    mrmse{ii} = median(feGetRep(eval(sprintf('%s',fes{ii})),'vox rmse'));
    mr2{ii}   = median(feGetRep(eval(sprintf('%s',fes{ii})),'vox r2'));
    dsd{ii}   = var(feGetRep(eval(sprintf('%s',fes{ii})),'dsigdemeaned'));
    dr2{ii}   = median(feGetRep(eval(sprintf('%s',fes{ii})),'vox r2 data'));
    mrp{ii}   = median(feGetRep(eval(sprintf('%s',fes{ii})),'voxelr2pearson'));
    
end

mrvNewGraphWin('Data RMSE (data reliability)')
for ii=1:length(os)
    bar(ii,drmse{ii},'FaceColor',colors{ii});
    hold on
end
set(gca,'ylim',[28 30],'ytick',[28 29 30], ...
    'xtick',[1 2 3],'xticklabel',[os],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('RMSE');xlabel('Algorithm type')

mrvNewGraphWin('Comparison on R^2 across tractography models')
for ii=1:length(os)
    h(ii) = semilogx(1:length(r2{ii}),r2{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii})    ;
    hold on

end
set(gca,'ylim',[20 40],'ytick',[20 30 40],'xtick',[0 10 100 1000],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('Percent variance explained');xlabel('iteration number')

mrvNewGraphWin('Comparison on RMSE across tractography models')
for ii=1:length(os)
    h(ii) = semilogx(1:length(rmse{ii}),rmse{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii})   ;
    hold on
end
set(gca,'ylim',[26 32],'ytick',[26 30 32],'xtick',[0 10 100 1000],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('RMSE');xlabel('iteration number')
legend(h,os,'box','off')

mrvNewGraphWin('Comparison on rRMSE across tractography models')
for ii=1:length(os)
    h(ii) = semilogx(1:length(rrmse{ii}),rrmse{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii})   ;
    hold on
end
set(gca,'ylim',[0.75 1.25],'ytick',[.75 1 1.25],'xtick',[0 10 100 1000],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('R_{rmse}');xlabel('iteration number')
legend(h,os,'box','off')

mrvNewGraphWin('Comparison on number of deleted fibers across tractography models')
for ii=1:length(os)
    h(ii) = loglog(.0001+(1:length(nFibs{ii})),nFibs{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii});
    hold on
end
set(gca,'xlim',[1 1000],'ylim',[1000 64000],'ytick',[1000 2000 4000 8000 16000 32000 64000], ...
    'xtick',[1 10 100 1000],'tickdir','out','box','off','FontSize',fontSiz)
ylabel('Number of fascicles in connectome');xlabel('iteration number')
legend(h,os,'box','off')

mrvNewGraphWin('Deleted fibers vs. rmse')
for ii=1:length(os)
    h(ii) = semilogx(nFibs{ii},rmse{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii}) ;
    hold on
end
set(gca,'ylim',[26 32],'ytick',[26 29 32],'xlim',[1000 64000],'xtick',[1000 2000 4000 8000 16000 32000 64000],'tickdir','out','box','off','FontSize',fontSiz)
xlabel('Number of fascicles in connectome');ylabel('rmse')
legend(h,os,'box','off')

mrvNewGraphWin('Deleted fibers vs. R^2')
for ii=1:length(os)
    h(ii) = semilogx(nFibs{ii},r2{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii})   ;
    hold on
end
set(gca,'ylim',[20 40],'ytick',[20 30 40],'xlim',[1000 64000],'xtick',[1000 2000 4000 8000 16000 32000 64000],'tickdir','out','box','off','FontSize',fontSiz)
xlabel('Number of fascicles in connectome');ylabel('R^2')
legend(h,os,'box','off')

mrvNewGraphWin('R_{rmse} vs. R^2')
for ii=1:length(os)
    h(ii) = semilogx(rrmse{ii},r2{ii},'ko-','color',colors{ii},'MarkerFaceColor',colors{ii})   ;
    hold on
end
set(gca,'ylim',[20 40],'ytick',[20 30 40],'xlim',[0.9 1.1],'xtick',[0.95 1 1.05],'tickdir','out','box','off','FontSize',fontSiz)
xlabel('R_{rmse}');ylabel('R^2')
legend(h,os,'box','off')
