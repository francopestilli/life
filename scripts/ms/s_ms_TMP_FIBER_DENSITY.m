% Load an FE structure
clx

% Medium connectome
%load /azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep3/fe_structures/0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax2__stream-ructures/0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax2__stream-500000_diffModAx100Rd0_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax2__stream-500000_diffModAx100Rd0_.mat% Good connectome

% Good connectome
load /azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep3/fe_structures/0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax4__m_prob-500000_diffModAx100Rd0_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax4__m_prob-500000_diffModAx100Rd0_.mat

% Bad connectome
%load /azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep3/fe_structures/0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_dwi_0005_0_tensor-500000_diffModAx100Rd0_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_dwi_0005_0_tensor-500000_diffModAx100Rd0_.mat

% get the fiber density
fd = feGet(fe,'fiber density');

% R rmse
R = feGetRep(fe,'vox rmse ratio');

% rmse
rmse = feGetRep(fe,'vox rmse');

% mean voxel signal
mdsig = feGetRep(fe,'meanvoxelsignal');

% non-diffusion weighted signal
s0 = feGetRep(fe,'b0vox');


% Make plots
mrvNewGraphWin('Fiber density');
set(gcf,'color','w')

% Fiber density before life
edges = logspace(.5,3.2,100);
centers = sqrt(edges(1:end-1).*edges(2:end));
y = histc(fd,edges)/length(fd)*100;
h = bar(y,'r');
set(h,'edgecolor','r','linewidth',.01)
set(get(h,'Children'),'FaceAlpha',.5,'EdgeAlpha',.5)
set(gca,'ylim',[0 3],'xlim',[1 100])
ticks = get(gca,'xtick');

title('Full connectome')
xlabel('Number of fibers per voxel')
set(gca,  'ytick',[0 1 2 3], ...
  'xticklabel', ceil(centers(ticks-1)) ,...
  'box','off','tickDir','out','xscale','lin')

% Fiber density after life
hold on
y = histc(fdw,edges)/length(fd)*100;
h = bar(y,'b');
set(h,'edgecolor','b','linewidth',.01)
set(get(h,'Children'),'FaceAlpha',.35,'EdgeAlpha',.35)
set(gca,'ylim',[0 3],'xlim',[1 100])
ticks = get(gca,'xtick');

ylabel('Percent white-matter volume')
xlabel('Number of fibers per voxel')
set(gca,  'ytick',[0 1 2 3], ...
  'xticklabel', ceil(centers(ticks-1)) ,...
  'box','off','tickDir','out','xscale','lin')

% Weigth density after life
mrvNewGraphWin('Weight density');
set(gcf,'color','w')

edges = logspace(-7,0,100);
centers = sqrt(edges(1:end-1).*edges(2:end));
y = histc(fw(:,2),edges)/length(fd)*100;
h = bar(y,'k');
set(h,'edgecolor','k','linewidth',.01)
title('Weights'' fiber density')
ylabel('Percent white-matter volume')
xlabel('Total fascicle''s contribution to the signal')
set(gca,'ylim',[0 6],'xlim',[.5 100])
ticks = get(gca,'xtick');

set(gca, 'ytick',[0 3 6], ...
  'xticklabel', ceil(1000000*centers(ticks-1))/1000000 ,...
  'box','off','tickDir','out','xscale','lin')

% Weigth density after life
mrvNewGraphWin('R_rmse');
set(gcf,'color','w')

edges = logspace(-.3,.6,100);
centers = sqrt(edges(1:end-1).*edges(2:end));
y = histc(R,edges)/length(R)*100;
h = bar(y,'k');
set(h,'edgecolor','k','linewidth',.01)
ylabel('Percent white-matter volume')
xlabel('R_{rmse}')
set(gca,'ylim',[0 6],'xlim',[.5 100])
ticks = get(gca,'xtick');

set(gca, 'ytick',[0 3 6], ...
  'xticklabel', ceil(100*centers(ticks-1))/100 ,...
  'box','off','tickDir','out','xscale','lin')

% Weigth density after life
mrvNewGraphWin(sprintf('R_{rmse} vs. FW scatter'));
set(gcf,'color','w')
loglog([nanmedian(fw(:,2)) nanmedian(fw(:,2))],[0.5 4],'k-')
hold on
plot([10^-8 10^0],[nanmedian(R) nanmedian(R)],'k-')
plot([10^-8 10^0],[1 1],'k--')
plot(fw(:,2),R,'ro','MarkerFaceColor','r')

xlabel('Sum of weights in voxel')
ylabel('R_{rmse}')
set(gca,'ylim',[0.5 4], 'ytick',[.5 1 2 4], ...
'box','off','tickDir','out')
axis square

% Weigth density after life
mrvNewGraphWin(sprintf('RMSE vs. FW scatter'));
set(gcf,'color','w')
semilogx([nanmedian(fw(:,2)) nanmedian(fw(:,2))],[0 120],'k-')
hold on
plot([10^-8 10^0],[nanmedian(rmse) nanmedian(rmse)],'k-')
plot(fw(:,2),rmse,'ro','MarkerFaceColor','r')

xlabel('Sum of weights in voxel')
ylabel('rmse')
set(gca,'ylim',[0 120], 'ytick',[0 60 120], ...
'box','off','tickDir','out')
axis square

% Weigth density after life
mrvNewGraphWin(sprintf('RMSE vs. R scatter'));
set(gcf,'color','w')
semilogx([nanmedian(R) nanmedian(R)],[0 120],'k-')
hold on
plot([10^-1 10^0],[nanmedian(rmse) nanmedian(rmse)],'k-')
plot(R,rmse,'ro','MarkerFaceColor','r')

xlabel('R_{rmse}')
ylabel('rmse')
set(gca,'ylim',[0 120], 'ytick',[0 60 120], ...
'box','off','tickDir','out')
axis square

% Weigth density after life
mrvNewGraphWin('R_rmse vs. FDW scatter');
set(gcf,'color','w')
loglog([nanmedian(fdw) nanmedian(fdw)],[0.5 4],'k-')
hold on
plot([10^-1 10^4],[nanmedian(R) nanmedian(R)],'k-')
plot([10^-1 10^4],[1 1],'k--')
plot(fdw,R,'ro','MarkerFaceColor','r')

xlabel('Fiber density after LiFE')
ylabel('R_{rmse}')
set(gca,'ylim',[0.5 4], 'ytick',[.5 1 2 4], ...
'box','off','tickDir','out')
axis square

% Weigth density after life
mrvNewGraphWin('R_rmse vs. FD scatter');
set(gcf,'color','w')
loglog([nanmedian(fd) nanmedian(fd)],[0.5 4],'k-')
hold on
plot([10^-1 10^4],[nanmedian(R) nanmedian(R)],'k-')
plot([10^-1 10^4],[1 1],'k--')
plot(fd,R,'ro','MarkerFaceColor','r')

xlabel('Fiber density before LiFE')
ylabel('R_{rmse}')
set(gca,'ylim',[0.5 4], 'ytick',[.5 1 2 4], ...
'box','off','tickDir','out')
axis square

% Mean dsig vs. sum fiber weights
mrvNewGraphWin('mean DSIG vs. FDW scatter');
set(gcf,'color','w')
loglog([nanmedian(fw(:,2)) nanmedian(fw(:,2))],[10^1 10^3],'k-')
hold on
plot([10^-7 10^0],[nanmedian(mdsig) nanmedian(mdsig)],'k-')
plot(fw(:,2),mdsig,'ro','MarkerFaceColor','r')

xlabel('Fiber density before LiFE')
ylabel('R_{rmse}')
set(gca,'ylim',[0.5 4], 'ytick',[.5 1 2 4], ...
'box','off','tickDir','out')
axis square

% Mean dsig vs. mean fiber weights
mrvNewGraphWin('mean DSIG vs. mean FW scatter');
set(gcf,'color','w')
loglog([nanmedian(fw(:,3)) nanmedian(fw(:,3))],[10^1 10^3],'k-')
hold on
plot([10^-7 10^0],[nanmedian(mdsig) nanmedian(mdsig)],'k-')
plot(fw(:,3),mdsig,'ro','MarkerFaceColor','r')
xlabel('Mean of fiber weight')
ylabel('Mean dSig')
set(gca, ...
'box','off','tickDir','out')
axis square

% Mean dsig vs. mean fiber weights
mrvNewGraphWin('DSIG vs. variance FW scatter');
set(gcf,'color','w')
loglog([nanmedian(fw(:,4)) nanmedian(fw(:,4))],[10^1 10^3],'k-')
hold on
plot([10^-17 10^0],[nanmedian(mdsig) nanmedian(mdsig)],'k-')
plot(fw(:,4),mdsig,'ro','MarkerFaceColor','r')
xlabel('Var of fiber weight')
ylabel('Mean dSig')
set(gca, 'box','off','tickDir','out')
axis square

% Mean dsig vs. mean fiber weights
mrvNewGraphWin('mean vs. variance FW scatter');
set(gcf,'color','w')
loglog([nanmedian(fw(:,4)) nanmedian(fw(:,4))],[10^-8 10^-1],'k-')
hold on
plot([10^-15 10^0],[nanmedian(fw(:,3)) nanmedian(fw(:,3))],'k-')
plot(fw(:,4),fw(:,3),'ro','MarkerFaceColor','r')
xlabel('Var of fiber weight')
ylabel('Mean of fiber weight')
set(gca, 'box','off','tickDir','out')
axis square

% Mean dsig vs. mean fiber weights
mrvNewGraphWin('mean dsig vs. s0 scatter');
set(gcf,'color','w')
loglog([nanmedian(mdsig) nanmedian(mdsig)],[10^2 10^3],'k-')
hold on
plot([10^2 10^3],[nanmedian(s0) nanmedian(s0)],'k-')
plot(mdsig,s0,'ro','MarkerFaceColor','r')
xlabel('Mean diffusion signal')
ylabel('S0 signal')
set(gca, 'box','off','tickDir','out')
axis square

% Mean dsig vs. mean fiber weights
mrvNewGraphWin('mean dsig vs. s0 scatter');
set(gcf,'color','w')
%loglog([nanmedian(mdsig) nanmedian(mdsig)],[10^2 10^3],'k-')
%hold on
%plot([10^2 10^3],[nanmedian(s0) nanmedian(s0)],'k-')
semilogx(fw(:,4),fw(:,2),'ro','MarkerFaceColor','r')
xlabel('Mean diffusion signal')
ylabel('S0 signal')
set(gca, 'box','off','tickDir','out')
axis square

% End
