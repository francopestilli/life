function fePlotQualityOfFit(fe)
% 
% function mctPlotQualityOfFit(r2_life,r2_orig,pSig,dSig)
%
% This function plots the results of a LiFE fit.
%
% It receives as inputs the R2 for the LiFe solution and for the original
% fiber model, the signal predicted by LiFE and the original signal.
%
% Eventually it will return a map of r2 for each voxel.
%
% Franco
%
% (C) 2012 Stanford VISTA team.

%% Do some feGets here

r2_life   = feGet(fe,'r2 life');
r2_orig   = feGet(fe,'r2 orig');
% pSig_life = feGet(fe,'psig life');
% pSig_orig = feGet(fe,'psig orig');
% dSig      = feGet(fe,'dsig');

%% Colors
cLiFE = [.4576 .5424 1];cOrig = [1 .5 0];

% plot the coefficient of determination
%subplot(1,2,1)
patch([.6 .6 1.5 1.5],[0 r2_life, r2_life, 0],'r','EdgeColor','w','FaceColor',cLiFE)
hold on
patch([1.6 1.6 2.5 2.5],[0 r2_orig r2_orig 0],'k','EdgeColor','w','FaceColor',cOrig)
text(.75,r2_life+3,sprintf('%2.2f',r2_life),'FontSize',16);
text(1.75,r2_orig+3,sprintf('%2.2f',r2_orig),'FontSize',16);

set(gca, ...
  'yTick',[0 25 50 75 100], ...
  'yLim',[-.1 101], ...
  'xTick',[.5 1 2 2.5], ...
  'xTickLabel',{'','LiFE','Original',''},...
  'xLim',[.5 3], 'TickDir','out')
ylabel('Percent variance explained (R^2)')
xlabel('White-matter fascicles.')

% Plot the predicted and observed signal
%subplot(1,2,2)
%lim = [min([dSig;pSig_orig]) max([dSig; pSig_orig])];
%plot(lim,lim,'k-'), hold on
%plot(dSig,pSig_orig,'k.','Color',cOrig);
%plot(dSig,pSig_life,'r.','Color',cLiFE);
%set(gca,'TickDir','out','Box','off','XLim',lim,'YLim',lim);
%axis('square');
%ylabel('Predicted diffusion signal')
%xlabel('Measured diffusion signal')

end
