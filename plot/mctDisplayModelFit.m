function h = mctDisplayModelFit(wFibers,wIso,dSig,life_pSig,life_r2, orig_pSig,orig_r2,fitType)
% 
% function h = mctDisplayModelFit(wFibers,wIso,fitType)
%
% Display the microtrack fits, still working on this. 
%
% Franco
%
% (C) 2012 Stanford VISTA team. 

h = mrvNewGraphWin(sprintf('LiFE Results - Using: %s',fitType));

% Fiber weights
cLiFE = [.4576 .5424 1];cOrig = [1 .5 0];
subplot(2,2,1), plot(wFibers,'o','MarkerFaceColor',cLiFE,'Color',cLiFE), axis tight
set(gca,'yLim',[-.1 max(wFibers) + 0.01],'TickDir','out','box','off');
title('LiFE solution - Fibers weights')
xlabel('Fiber number');ylabel('Fiber weight')

% isotropic components' weights, mean signal in a voxel
subplot(2,2,2), plot(wIso,'k.'), axis tight
set(gca,'YLim',[-0.1 nanmax(wIso(:))*1.1],'TickDir','out','box','off');
title('Mean DW signal in each voxel')
xlabel('Voxel number / diffusion direction index');ylabel('DW signal')

% Plot the signals in the voxels
subplot(2,2,3:4) 
plot(dSig,'k-'); hold on
plot(life_pSig,'-','LineWidth',1,'MarkerFaceColor',cLiFE,'Color',cLiFE,'MarkerEdgeColor','w');
plot(orig_pSig,'--','LineWidth',1,'MarkerFaceColor',cOrig,'Color',cOrig,'MarkerEdgeColor','w');
axis tight
xlabel('Voxels and diffusion directions');
ylabel('DW signal')
title(sprintf('Percent variance explained: Orig:%2.2f, LiFE:%2.2f',orig_r2,life_r2))
set(gca,'TickDir','out','box','off');

return
