function s_life_tensor(savemovie,moviedir)
%
% Script used to visualize the tensor as a 3D and 3D signal.
% Saves movies
%
% Franco
%
% Chosen LiFE colors
% cLiFE = [.4576 .5424 1];cOrig = [1 .5 0];
global gPlotBox

if notDefined('moviedir'),moviedir   = '/biac2/wandell6/data/frk/LiFE/results/mpi/tensor_movies';end
if notDefined('figdir'),figdir   = '/biac2/wandell6/data/frk/LiFE/results/mpi/tensor_figs';end

if notDefined('savemovie'),savemovie = 0;end
if notDefined('whichTensor'),whichTensor = 'long';end

gPlotBox     = 0;
plotDirs     = 0;
doPlotADC    = 0;
doPlotDSig   = 1;
doPlotTensor = 0;
doPlotSigOnTensor = 0;

% Pick a voxel
switch whichTensor
  case {'s','spherical',}
    %coords = [47 54 43];  % Circular
    coords = [55 76 38];  % Circular
  case {'cigar','long','ellipsoidal'}
    %coords = [44 54 43];  % Directional
    coords = [41 43 38];  % Directional
  otherwise
    keyboard
end

%% Load diffusion weighted imaging data
% The vistadata diffusion sample data are 40-directions.  The directory
% contains the dwi data as well as the bvals and bvecs. 
dataDir = fullfile(mrvDataRootPath,'diffusion','sampleData');
dwi = dwiLoad(fullfile(dataDir,'raw','dwi.nii.gz'));

% Plot the diffusion directions measured in the scanner (bvecs) on a
% sphere.
if plotDirs
  [fh ah] = plotDirections(dwi);
  if savemovie
    moviename = fullfile(moviedir,sprintf('%s','bvecs'));
    makeMovie(fh,ah,moviename);
  end
end


%% Compute and Plot the ADC, Apparent Diffusion Coefficients
% Compute the diffusion-weighted signal and the apparent diffusion
% coefficients
[dSig ADC] =  computeDSigADC(dwi,coords);
%
% Equivalent to
%
% dSig = dwiGet(dwi,'diffusion data image',coords);
% ADC = dwiGet(dwi,'adc data image',coords);
%

% compute the quadratic form of the tensor
Q = dwiTensor(dwi,ADC);

if doPlotADC
  % Plot the apparent diffusion coefficiets - 3D
  [fh ah] = plotADC(dwi,ADC,Q);
  
  if savemovie
    moviename = fullfile(moviedir,sprintf('%s','ADC_3D_peanut'));
    makeMovie(fh,ah,moviename);
  end
end

% Now use the tensor model to predict the ADC given the measured ADC
predADC = predictedADC(dwi,ADC,Q);

% Plot the predicted vs. the measured ADC signal.
fh = plotMeasuredPredictedADC(predADC,ADC,figdir);

% Plot the series of ADC
fh = plotADCseries(ADC,predADC,figdir);

%% Plot Tensor, Signal measured and predicted 3D, 3D and as a series
% We have a prediction ADC = u'Qu.  The diffusion ellipsoid associated with
% Q are the scaled unit vectors such that v'inv(Q)v = 1.  The vectors v =
% sqrt(ADC)*u will be close to the ellipsoid defined by Q.

if doPlotDSig
  % Plot Plot Signal and best Fit tensor model in 3D
  [fh ah] = plotDSig(dwi,dSig,coords,Q);
  
  if savemovie
    moviename = fullfile(moviedir,sprintf('%s','DSIG_3D_squash'));
    makeMovie(fh,ah,moviename);
  end
end

if doPlotTensor
  % Plot Tensor with predictd and measured signal in 3D
  [fh ah] = plotTensorWithDSig(dwi,dSig,coords,Q,doPlotSigOnTensor);
  
  if savemovie
    moviename = fullfile(moviedir,sprintf('%s','TENSOR_3D_cigar'));
    makeMovie(fh,ah,moviename);
  end
end

% Plot compute the predicted signal and plot it in 2D
predDSig = predictedDSig(dwi,coords,Q);
h = plotMeasuredPredictedDSig(predDSig,dSig,figdir);

% Plot the series of ADC
h = plotDSigSeries(dSig,predDSig,figdir);

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A comparison of the estimated and observed ADC values
mrvNewGraphWin;
subplot(1,2,1)
err = 100*(predADC - ADC) ./ ADC;
hist(err)
xlabel('Percent Error')
ylabel('N voxels')

mean(err)
std(err)

err3 = diag(err)*bvecs;
subplot(1,2,2)
plot3(err3(:,1),err3(:,2),err3(:,3),'.')
axis equal

%% A comparison of the estimated and observed signal values

%
%    d = S0 exp(-b (u'Qu))
%
dSigEst = S0*exp(-bvals.* diag(bvecs*Q*bvecs'));

mrvNewGraphWin
subplot(1,2,1)
plot(dSigEst,dSig,'.')
xlabel('Estimated diffusion sig')
ylabel('Measured diffusion sig')

subplot(1,2,2), 
p = 100*(dSigEst(:) - dSig(:)) ./ dSig(:);
hist(p,20);
xlabel('% Error')

%% The error in 3-space
p3 = diag(p)*bvecs;
plot3(p3(:,1),p3(:,2),p3(:,3),'.')
axis equal
axis square
set(gca,'TickDir','out','XTick',[],'YTick',[],'ZTick',[])

end

%%%%%%%%%%%%%%%%%%
% plotDSigSeries %
%%%%%%%%%%%%%%%%%%
function fh = plotDSigSeries(dSig,predDSig,moviedir)
fh = mrvNewGraphWin('Diffusion signal series');
whitebg(fh,[1 1 1])

col1 = [.4576 .5424 1];
col2 = [1 .5 0];

plot(predDSig,'w-','Color',col2,'LineWidth',2)
hold on
plot(dSig,'co','MarkerSize',12,'MarkerFaceColor',[.1 .1 .1],'MarkerEdgeColor','w')
r2 = calccod(dSig(:),predDSig(:));
r2 = 100*corr(dSig(:),predDSig(:))^2;
title(sprintf('Variance explained %2.2f%%',r2))
xlabel('Measured diffusion directions (list)')
ylabel('Diffusion signal')
set(gca,'TickDir','out','Box','off');
set(fh,'Color',[1 1 1]);
box off  
set(gcf,'Position',[485 848 1397 250])

%savefigvista(fh ,'diffusion_signal_series','eps',moviedir,1,0);

end



%%%%%%%%%%%%%%%%%
% plotADCseries %
%%%%%%%%%%%%%%%%%
function fh = plotADCseries(ADC,predADC,moviedir)
fh = mrvNewGraphWin('ADC series');
whitebg(fh,[1 1 1])

col1 = [.4576 .5424 1];
col2 = [1 .5 0];
plot(predADC,'w-','Color',col2,'LineWidth',2)
hold on
plot(ADC,'co','MarkerSize',14,'MarkerFaceColor',[.1 .1 .1],'MarkerEdgeColor','w')
r2 = calccod(ADC,predADC);
r2 = 100*corr(ADC(:),predADC(:))^2;
set(gca,'TickDir','out','Box','off');
set(fh,'Color',[1 1 1]);
set(gcf,'Position',[485 848 1397 250])

title(sprintf('Variance explained %2.2f%%',r2))
xlabel('Measured diffusion directions (list)')
ylabel('ADC (diffusion weighted)')
box off
%savefigvista(fh ,'adc_series','eps',moviedir,1,0);

end


%%%%%%%%%%%%%%%%%
% predictedDSig %
%%%%%%%%%%%%%%%%%
function predDSig = predictedDSig(dwi,coords,Q)
% Predict the diffusion signal.
S0 = dwiGet(dwi, 'S0 image',coords);
predDSig = dwiComputeSignal(S0, dwiGet(dwi,'diffusion bvecs'), dwiGet(dwi,'diffusion bvals'), Q(:)');
end


%%%%%%%%%%%%%%%%%%%%%%
% plotTensorWithDSig %
%%%%%%%%%%%%%%%%%%%%%%
function [fh ah] = plotTensorWithDSig(dwi,dSig,coords,Q,doPlotSigOnTensor)
global gPlotBox

fh = mrvNewGraphWin('Tensor and Diffusion Signal, 3D');
ellipsoidFromTensorLocal(Q);
whitebg([0 0 0])
hold on, grid off
set(fh,'Color',[0 0 0]);

if doPlotSigOnTensor
  % Predict the diffusion signal.
  S0 = dwiGet(dwi, 'S0 image',coords);
  predDSig = dwiComputeSignal(S0, dwiGet(dwi,'diffusion bvecs'), dwiGet(dwi,'diffusion bvals'), Q(:)');
  
  % Compute and plot vector of measured adcs
  if coords(1) == 41
    ds = diag(dSig)*dwiGet(dwi,'diffusion bvecs')/dwiGet(dwi,'b0image',coords) .*0.95;
    pds = diag(predDSig)*dwiGet(dwi,'diffusion bvecs')/dwiGet(dwi,'b0image',coords) .* 0.95;
    
  elseif coords(1) == 55
    ds = diag(dSig)*dwiGet(dwi,'diffusion bvecs')/dwiGet(dwi,'b0image',coords) .* 1.5;
    pds = diag(predDSig)*dwiGet(dwi,'diffusion bvecs')/dwiGet(dwi,'b0image',coords) .* 1.5;
    
  elseif coords(1) == 44
    ds = diag(dSig)*dwiGet(dwi,'diffusion bvecs')/dwiGet(dwi,'b0image',coords) .* 2;
    pds = diag(predDSig)*dwiGet(dwi,'diffusion bvecs')/dwiGet(dwi,'b0image',coords) .* 2;
    
  else
    ds = diag(dSig)*dwiGet(dwi,'diffusion bvecs')/dwiGet(dwi,'b0image',coords) .* 16;
    pds = diag(predDSig)*dwiGet(dwi,'diffusion bvecs')/dwiGet(dwi,'b0image',coords) .* 16;
    
    
  end
  
  plot3(ds(:,1),ds(:,2),ds(:,3),'o','MarkerSize',16,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.35 .35 .35])
  %plot3(pds(:,1),pds(:,2),pds(:,3),'d','MarkerSize',16,'MarkerFaceColor',[.9 .7 .4],'MarkerEdgeColor',[.35 .35 .35])
end

%title('Diffusion distance ellipsoid')
set(gca,'TickDir','out')
set(gca,'TickDir','out','XTick',[],'YTick',[],'ZTick',[])

if gPlotBox
  box on
else
  axis off
end

ah = gca;
end

  
%%%%%%%%%%%%
% plotDSig %
%%%%%%%%%%%%
function [fh ah] = plotDSig(dwi,dSig,coords,Q)
global gPlotBox

% Start the figure
fh = mrvNewGraphWin('DWI signal - 3D');
whitebg([0 0 0]), box on, grid off
set(fh,'Color',[0 0 0]);

grid off
% Plot the predicted signal ans a surface
if ~notDefined('Q')
  cmap = cool(255);
  [X,Y,Z] = sphere(60);
  [r,c] = size(X);  
  v        = [X(:),Y(:),Z(:)]; % simulatd bvecs for the surface
  bv       = unique(dwiGet(dwi,'bvals'));
  predDSig = dwiComputeSignal(1, v, bv(end).*ones(size(v,1),1), Q(:)');
  pds      = diag(predDSig)*v;
  
  v = pds;
  x = reshape(v(:,1),r,c);
  y = reshape(v(:,2),r,c);
  z = reshape(v(:,3),r,c);
  surf(x,y,z,repmat(256,r,c),'EdgeAlpha',0.1);
  axis equal, axis on, box on, grid off
  colormap([cmap; .25 .25 .25]), alpha(0.5)
  camlight; lighting phong; material shiny;
  set(gca, 'Projection', 'perspective');
  hold on
end

% Plot the Measured Signal
S0 = dwiGet(dwi, 'S0 image',coords);
ds = diag(dSig)*dwiGet(dwi,'diffusion bvecs')/dwiGet(dwi,'b0image',coords);
plot3(ds(:,1),ds(:,2),ds(:,3),'o','MarkerSize',12,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.35 .35 .35])
set(gca,'TickDir','out')
set(gca,'TickDir','out','XTick',[],'YTick',[],'ZTick',[])

if gPlotBox
  box on
else
  axis off
end
ah = gca;
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% plotMeasuredPredicted %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [fh ah] = plotMeasuredPredictedADC(ADCpred,ADC,moviedir)
% We would like to have a measure of how well the model does compared to a
% repeated measure.  Or, we would like a bootstrap of the estimate.  That
% is done in dtiRawFitTensor

fh = mrvNewGraphWin('Apparent Diffusion Coefficients - Predicted vs. Measured');
col1 = [.4576 .5424 1];
col2 = [1 .5 0];
plot([.2 4],[.2 4],'w-')
hold on
plot(ADC,ADCpred,'o','MarkerSize',12,'MarkerFaceColor',col2,'MarkerEdgeColor','w');
axis equal; grid off;
r2 = calccod(ADC,ADCpred);
r2 = 100*corr(ADC(:),ADCpred(:))^2;

title(sprintf('Variance explained %2.2f%%',r2))
xlabel('Estimated Apparent Diffusion Coefficients');
ylabel('Predicted Apparent Diffusion Coefficients');
set(gca,'TickDir','out','ylim',[.2 4],'xlim',[.2 4],'yscale','log','xscale','log', ...
         'ytick',[.2 .4 .8 1.6 3.2],'xtick',[.2 .4 .8 1.6 3.2])
%axis equal
axis square
whitebg([1 1 1])
set(gca,'TickDir','out','Box','off');
set(fh,'Color',[1 1 1]);
ah = gca;
%savefigvista(fh ,'scatter_ADC','eps',moviedir,1,0);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotMeasuredPredicteDSig %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fh ah] = plotMeasuredPredictedDSig(predDSig,dSig,moviedir)
% We would like to have a measure of how well the model does compared to a
% repeated measure.  Or, we would like a bootstrap of the estimate.  That
% is done in dtiRawFitTensor

fh = mrvNewGraphWin('MR Diffusion Signal - Predicted vs. Measured');
col1 = [.4576 .5424 1];
col2 = [1 .5 0];
hold on
plot([0 700],[0 700],'k-')
plot(dSig,predDSig,'o','MarkerSize',12,'MarkerFaceColor',col2,'MarkerEdgeColor','w');
r2 = calccod(dSig(:),predDSig(:));
r2 = 100*corr(dSig(:),predDSig(:))^2;

title(sprintf('Variance explained %2.2f%%',r2))
xlabel('Measured diffusion signal');
ylabel('Predicted diffusion signal');
set(gca,'TickDir','out','ylim',[0 700],'xlim',[0 700])
%axis equal
axis square
grid off
whitebg([1 1 1])
set(gca,'TickDir','out','Box','off');
set(fh,'Color',[1 1 1]);
ah = gca;
%savefigvista(fh ,'scatter_Sig','eps',moviedir,1,0);

end


%%%%%%%%%%%
% predADC %
%%%%%%%%%%%
function predADC = predictedADC(dwi,ADC,Q)
% Solve for the quadratic, Q, that predicts the ADC values. 
%
% The ADC values are measured.  Each measurement has a direction and
% amplitude, m = bvecs*bval.  We want to find Q such that 
%
%    (bvec(:)'* Q * bvec(:) = ADC
%
% We express this as a quadratic equation.  Suppose the entries of Q are
% (sorry about this) qij.  Suppose that for each direction, bvec, we have a
% particular level of b.  Then we have for each direction,
%
%  ADC(:) = q11*b(:,1)^2 + ... 2*qij*b(:,i)*b(:,j) + ... q33*b(:,3)^2
%
% The coefficients qij are the same in every equation.  So, we can pull
% them out into a column vector.  We have a matrix, V, whose entries are
% these b values.
%
%   ADC = V*q
%
bvecs = dwiGet(dwi,'diffusion bvecs');

% To compare the observed and predicted, do this
predADC = zeros(size(ADC));
for ii = 1:size(bvecs,1)
    u = bvecs(ii,:);
    predADC(ii) = u(:)'*Q*u(:);
end
end


%%%%%%%%%%%%%
% dwiTensor %
%%%%%%%%%%%%%
function Q = dwiTensor(dwi,ADC)

% We compute using only the diffusion-weighted data.
%
% Because the bvecsP enter the equation on both sides of Q, we take the
% square root of the bvals and multiply it times the bvecs.
% b = bvecs .* repmat(bvals.^0.5,1,3);  % Scaled, non-zero bvecs
% The diffusion weighted bvecs
bvecs = dwiGet(dwi,'diffusion bvecs');
b = bvecs;  % Do not scale.  We scale when we compute diffusion. 

% Here is the big matrix
V = [b(:,1).^2, b(:,2).^2, b(:,3).^2, 2* b(:,1).*b(:,2), 2* b(:,1).*b(:,3), 2*b(:,2).*b(:,3)];

% Now, we divide the matrix V by the measured ADC values to obtain the qij
% values in the parameter, tensor
tensor = V\ADC;

% We convert the format from a vector to a 3x3 Quadratic
Q = dt6VECtoMAT(tensor);  % eigs(Q)
% svd(Q)

end % of function here 


%%%%%%%%%%%
% plotADC %
%%%%%%%%%%%
function [fh ah] = plotADC(dwi,adc,Q)
global gPlotBox

% Start the figure
fh = mrvNewGraphWin('Apparent Diffusion Coefficients - 3D');
whitebg([0 0 0])
set(fh,'Color',[0 0 0]);

% We use this to get the predicted ADC values from the tensor
% We plot the surface if available.
if ~notDefined('Q')
  cmap = cool(255);

  % User passed in Q, make the predicted peanut
  [X,Y,Z] = sphere(60);
  [r,c] = size(X);
  
  v = [X(:),Y(:),Z(:)];
  adcPredicted = diag(v*Q*v');
  v = diag(adcPredicted)*v;
  
  x = reshape(v(:,1),r,c);
  y = reshape(v(:,2),r,c);
  z = reshape(v(:,3),r,c);
  surf(x,y,z,repmat(256,r,c),'EdgeAlpha',0.1);
  axis equal,
  colormap([cmap; .25 .25 .25]), alpha(0.5)
  camlight; 
  lighting phong; material shiny;
  %set(gca, 'Projection', 'perspective');
  hold on
end

% The diffusion weighted bvecs
bvecs = dwiGet(dwi,'diffusion bvecs');

% Compute and plot vector of measured adcs
adcV = diag(adc)*bvecs;
plot3(adcV(:,1),adcV(:,2),adcV(:,3),'o','MarkerSize',12,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.35 .35 .35])
grid off
set(gca,'TickDir','out')
set(gca,'TickDir','out','XTick',[],'YTick',[],'ZTick',[])

if gPlotBox
  box on
else
  axis off
end

ah = gca;
end


%%%%%%%%%%%%%%%%%%
% computeDSigADC %
%%%%%%%%%%%%%%%%%%
function [dSig ADC] =  computeDSigADC(dwi,coords)
% Find the non-difussion bvecs and bvals (b ~= 0).
bvals = dwiGet(dwi,'diffusion bvals');

% get the signal at B0 
S0 = dwiGet(dwi,'b0 image',coords);

% Diffusion data are nCoords x nImages, excludes the b=0 measures
dSig = dwiGet(dwi,'diffusion data image', coords);

% Calculate the ADC from the diffusion-weighted data
% dSig = S0 * exp(-b * ADC)
% ADC = bVec*Q*bVec'
%
ADC = - diag( (bvals).^-1 )*log(dSig(:)/S0);  % um2/ms
end


%%%%%%%%%%%%%
% makeMovie %
%%%%%%%%%%%%%
function m = makeMovie(fh,ah,filename)

% rotate the sphere and collect frames.
c = 1;
%rotations = [10 12 14 16 20 25 30:10:70 75 80 84 86 88 90];
%rotations = [rotations, fliplr(rotations)];
rotations = (-13:.4:13).^2 + 1;
   
set(fh,'Position',[0 0 1200 1600],'Visible','off');
axis off;
axis vis3d
set(gcf,'Color',[0 0 0]);
%set(gca,'LooseInset',get(gca,'TightInset'))
    
% Open a video object
vidObj = VideoWriter(filename);
vidObj.FrameRate = 15;
vidObj.Quality = 100;
open(vidObj);

for ii = rotations 
   %view(ah,ii,10); 
   
   camorbit(0.9,-0.1)
   drawnow
   m(c) = getframe(fh);
   writeVideo(vidObj,m(c));
   c = c + 1;
end
close(vidObj);
  
end


%%%%%%%%%%%%%%%%%%
% plotDirections %
%%%%%%%%%%%%%%%%%%
function [fh ah] = plotDirections(dwi,visible)
% To become a dwiGet(dwi,'bvals positive')
% Find the positive values

bvals = dwi.bvals;
bvecs = dwi.bvecs;
lst   = (bvals == 0);
bvals = bvals(~lst);
bvecs = bvecs(~lst,:);

tmp = diag(bvals)*bvecs; tmp = unique(tmp,'rows');
T = DelaunayTri(tmp(:,1),tmp(:,2),tmp(:,3));

fh = mrvNewGraphWin('Measured diffusion directions');
tetramesh(T,'Marker','o', 'MarkerFaceColor',[.8 .6 .2], ...
     'FaceAlpha',0.35,'FaceLighting','gouraud', ...
     'EdgeAlpha',0.001,'EdgeColor',[.5 .5 .5],'EdgeLighting','phong', ...
     'SpecularColorReflectance',.5, 'MarkerSize',12);
axis equal; axis off
set(gca,'Color',[0 0 0])
set (gcf,'Color',[0 0 0])
colormap(bone)
set(gca,'TickDir','out')
set(gca,'TickDir','out','XTick',[],'YTick',[],'ZTick',[])
ah = gca;

end


%%%%%%%%%%%%%%%%%%%%%%%
% ellipsoidFromTensor %
%%%%%%%%%%%%%%%%%%%%%%%
function ellipsoidFromTensorLocal(Q,C)
%Plot a diffusion ellipsoid from a tensor
%
% This is a copy of ellipsoid for tensor for local edits.
%
% (c) Stanford VISTA Team
if notDefined('N'), N = 60; end      % number of points on the sphere
if notDefined('C'), C = [0,0,0]; end % coordinates of the points in 3D space

% Make unit vectors on a sphere, u
[x,y,z] = sphere(N);   % build a sphere sampled with N resolution
u = [x(:), y(:), z(:)];% reorganize the univectors on the sheres (basically the diffusion directions scaled to the unit sphere)

sphAdc = zeros(size(u,1),1);
for ii=1:size(u,1)
  % Notice that this is multiplying by inv(Q), not Q.
  % Uncomment to see the equivalence of
  % sphADC(ii) = u(ii,:)*inv(Q)*u(ii,:)';
  sphAdc(ii) = u(ii,:)*(Q\u(ii,:)');
end

% The sqrt ADC is multipled into the unit vector because the factor is
% multiplied in twice (left and right side of inv(Q)).
v = diag(1 ./ sphAdc.^0.5)*u;  %Unit vectors in the rows.  Each row scaled.

% Now reshape and plot
sz = size(x);
x = reshape(v(:,1),sz);
y = reshape(v(:,2),sz);
z = reshape(v(:,3),sz);
x = x + C(1); y = y + C(2); z = z + C(3);

% No return arguments, so make the plot.
cmap = cool(255);

surf(x,y,z,repmat(256,size(z)),'EdgeAlpha',0.1);

axis equal, colormap([cmap; .25 .25 .25]), alpha(0.5)
camlight; lighting phong; material shiny;
set(gca, 'Projection', 'perspective');
set(gca,'TickDir','out')
set(gca,'TickDir','out','XTick',[],'YTick',[],'ZTick',[])


end

