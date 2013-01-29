
function [uData, g] = fePlot(fe,plotType,varargin)
% Gateway routine for plotting from the fascicle evaluation (fe) struct
%
%   [uData, g] = fePlot(fe,plotType,varargin)
%
% Example
%   fePlot(fe,'quality of fit');
%   fePlot(fe,'A matrix');
%
%   fePlot(fe,'bvecs');
%
%   fePlot(fe,'dsig measured');
%
% Franco (c) Stanford VISTA Team, 2012
%
%-----------
% Plot the fiber group in a matlab 3D mesh.
% uData = fePlot(fe,'connectome')
%-----------
% Plots the quality of fit of the LiFE model
% fData = ePlot(fe,'quality of fit');
%-----------
% Show an image of the LiFE model/matrix
% uData = fePlot(fe,'model')
%-----------
% Plot the series of the measured diffusion signal reordered in such a way
% that close angle in on the sphere are close together in the series. This
% call uses the functionality in feGet(fe,dsig measured', varargin)
% uData = fePlot(fe,'dsig measured')
% uData = fePlot(fe,'dsig measured',voxelIndex)
% uData = fePlot(fe,'dsig measured',coords)
%-----------
% Shows the bvecs on a sphere.
% uData = dwiPlot(fe,'bvecs')
%-----------
% Plot a histogram of the length of the fibers in the connectome.
% uData = fePlot(fe,'fiber length')
%-----------
%
%-----------
if notDefined('fe'),       error('''fe'' structure required.'); end
if notDefined('plotType'), error('plotType required'); end

% Open up your window for the current figure
g = mrvNewGraphWin(sprintf('LiFE - %s - %s',upper(plotType),feGet(fe,'name')));
uData = [];

% Format the input parameter
plotType = mrvParamFormat(plotType);

% Find the plot
switch plotType
  case {'map','maprep'}
    % Generate the requested map and save it to volume.
    if (length(plotType) > 3)
      % Use the xvalidation 'get' routine
      map = feValues2volume(feGetRep(fe, varargin{1}), ...
        feGet(fe,'roi coords'), ...
        feGetRep(fe,'map size'));
    else
      % use the default 'get' routine.
      map = feValues2volume(feGet(fe, varargin{1}), ...
        feGet(fe,'roi coords'), ...
        feGet(fe,'map size'));
    end
    
    % Get the slices to display (acpc) -120 to 120
    if (length(varargin)==1), slices = [-11:13];
    else slices = varargin{2};end
    
    % The xform is taken from the LiFE structure.
    mapXform = feGet(fe,'xform img2acpc');
    
    % Get the T1 anatomy from the session.
    ni         = niftiRead(feGet(fe,'anatomy file'));
    anatomyImg = double(ni.data);
    
    % Get the anatomy xform, from the file?
    anatomyXform = ni.qto_xyz;
    clear ni;
    
    % Get the range of the map to display from the data.
    clippingRange = [0 max(map(:))];%[max([0,min(map(:))]) max([max(map(:)),1])];
    mrAnatOverlayMontage(map, mapXform, anatomyImg, anatomyXform, [], clippingRange,slices);
    colormap('autumn');colorbar('location','eastoutside')
    title(upper(varargin{1}))
    
  case 'rmse'
    % Plot a histogram of the RMSE of the LiFE fit.
    % uData = fePlot(fe,'rmse')    
    uData.totVarExp = feGet(fe,'totpve');
    uData.totrmse   = feGet(fe,'totalrmse');
    uData.rmsebyVox   = feGet(fe,'voxrmse');
    uData.rmsemedian = nanmedian(uData.rmsebyVox);

    % male a histogram of the length of the fibers
    [y,x] = hist(uData.rmsebyVox,ceil(length(uData.rmsebyVox)/10));
    bar(x,y,.65,'r')
    hold on
    plot([uData.rmsemedian uData.rmsemedian],[0 400],'k')
    title(sprintf('Total percent variance explained: %2.2f\nTotal RMSE: %2.2f',uData.totVarExp, uData.totrmse))
    xlabel('Root mean square error')
    ylabel('Number of occurrences')
     
  case 'pve'
    % Plot a histogram of the percent variance explained in each voxel with
    % the LiFE fit.
    % uData = fePlot(fe,'pve')    
    uData.totVarExp = feGet(fe,'totpve');
    uData.totrmse   = feGet(fe,'totalrmse');
    uData.pve   = 100*feGet(fe,'voxr2zero');
    uData.pve( uData.pve < -1000 ) = -200; % if some of the values are -Inf set them to -200
    uData.medianpve = nanmedian(uData.pve);

    % male a histogram of the length of the fibers
    [y,x] = hist(uData.pve,ceil(length(uData.pve)/10));
    bar(x,y,1,'r')    
    hold on
    plot([uData.medianpve uData.medianpve],[0 400],'k')
    title(sprintf('Percent variance explained: %2.2f\nRMSE: %2.2f',uData.totVarExp, uData.totrmse))
    xlabel('Percent variance explained')
    ylabel('Number of occurrences')

  case 'rmserep'
    % Plot a histogram of the RMSE of the LiFE fit.
    % uData = fePlot(fe,'rmse')    
    uData.totVarExp  = feGetRep(fe,'totpve');
    uData.totrmse    = feGetRep(fe,'totalrmse');
    uData.rmsebyVox  = feGetRep(fe,'voxrmse');
    uData.rmsemedian = nanmedian(uData.rmsebyVox);
    
    % male a histogram of the length of the fibers
    [y,x] = hist(uData.rmsebyVox,ceil(length(uData.rmsebyVox)/10));
    bar(x,y,.65,'r')
    hold on
    plot([uData.rmsemedian uData.rmsemedian],[0 400],'k')
    title(sprintf('Total percent variance explained: %2.2f\nTotal RMSE: %2.2f',uData.totVarExp, uData.totrmse))
    xlabel('Root mean square error')
    ylabel('Number of occurrences')
     
  case 'pverep'
    % Plot a histogram of the percent variance explained in each voxel with
    % the LiFE fit.
    % uData = fePlot(fe,'pve')    
    uData.totVarExp = feGetRep(fe,'totpve');
    uData.totrmse   = feGetRep(fe,'totalrmse');
    uData.pve       = 100*feGetRep(fe,'voxr2zero');
    uData.pve( uData.pve < -1000 ) = -200; % if some of the values are -Inf set them to -200
    uData.medianpve = nanmedian(uData.pve);

    % male a histogram of the length of the fibers
    [y,x] = hist(uData.pve,ceil(length(uData.pve)/10));
    bar(x,y,1,'r')    
    hold on
    plot([uData.medianpve uData.medianpve],[0 400],'k')
    title(sprintf('Percent variance explained: %2.2f\nRMSE: %2.2f',uData.totVarExp, uData.totrmse))
    xlabel('Percent variance explained')
    ylabel('Number of occurrences')
     
  case 'rmserepdata'
    % Plot a histogram of the RMSE of the LiFE fit.
    % uData = fePlot(fe,'rmse')    
    uData.totVarExp  = feGetRep(fe,'totpvedata');
    uData.totrmse    = feGetRep(fe,'totalrmsedata');
    uData.rmsebyVox  = feGetRep(fe,'voxrmsedata');
    uData.rmsemedian = nanmedian(uData.rmsebyVox);
     
    % male a histogram of the length of the fibers
    [y,x] = hist(uData.rmsebyVox,ceil(length(uData.rmsebyVox)/10));
    bar(x,y,.65,'g')
    hold on
    plot([uData.rmsemedian uData.rmsemedian],[0 400],'k')
    title(sprintf('Total percent variance explained by the data : %2.2f\nTotal RMSE: %2.2f',uData.totVarExp, uData.totrmse))
    xlabel('Root mean square error')
    ylabel('Number of occurrences')
  
  case 'pverepdata'
    % Plot a histogram of the percent variance explained in each voxel with
    % the LiFE fit.
    % uData = fePlot(fe,'pve')    
    uData.totVarExp = feGetRep(fe,'totpvedata');
    uData.totrmse   = feGetRep(fe,'totalrmsedata');
    uData.pve       = 100*feGetRep(fe,'voxr2zerodata');
    uData.pve( uData.pve < -1000 ) = -200; % if some of the values are -Inf set them to -200
    uData.medianpve = nanmedian(uData.pve);
    
    % male a histogram of the length of the fibers
    [y,x] = hist(uData.pve,ceil(length(uData.pve)/10));
    bar(x,y,1,'g')    
    hold on
    plot([uData.medianpve uData.medianpve],[0 400],'k')
    title(sprintf('Percent variance explained by the data: %2.2f\nRMSE: %2.2f',uData.totVarExp, uData.totrmse))
    xlabel('Percent variance explained')
    ylabel('Number of occurrences')
  
  case 'rmseratio'
    % Plot a histogram of the RMSE of the LiFE fit.
    % uData = fePlot(fe,'rmseratio')
    uData.totVarExpD   = feGetRep(fe,'totpvedata');
    uData.totVarExpM  = feGetRep(fe,'totpve');
    uData.totVarExpMvw= 100*feGetRep(fe,'totalr2voxelwise');
    uData.totrmse     = feGetRep(fe,'totalrmseratio');
    uData.rmsebyVox   = feGetRep(fe,'voxrmseratio');
    uData.medianRmse  = nanmedian(uData.rmsebyVox);
    
    % male a histogram of the length of the fibers
    [y,x] = hist(uData.rmsebyVox,ceil(length(uData.rmsebyVox)/10));
    bar(x,y,.65,'k')
    hold on
    plot([uData.medianRmse uData.medianRmse],[0 400],'k')
    title(sprintf('PVE Data1/Data2 (%2.2f) | Global (%2.2f) | Voxel-wise (%2.2f)\nRMSE Ratio, Total (%2.2f) | Median (%2.2f)', ...
      uData.totVarExpD, ...
      uData.totVarExpM, ...
      uData.totVarExpMvw, ...
      uData.totrmse,   ...
      uData.medianRmse))
    xlabel('Root mean square error')
    ylabel('Number of occurrences')

  case 'rmseratiovoxelwise'
    % Plot a histogram of the RMSE of the LiFE fit.
    % uData = fePlot(fe,'rmseratio')
    uData.totVarExpD  = feGetRep(fe,'totpvedata');
    uData.totVarExpM  = feGetRep(fe,'totpve');
    uData.totVarExpMvw= 100*feGetRep(fe,'totalr2voxelwise');
    uData.totrmse     = feGetRep(fe,'totalrmseratiovoxelwise');
    uData.rmsebyVox   = feGetRep(fe,'voxrmseratiovoxelwise');
    uData.medianRmse  = nanmedian(uData.rmsebyVox);
    
    % male a histogram of the length of the fibers
    [y,x] = hist(uData.rmsebyVox,ceil(length(uData.rmsebyVox)/10));
    bar(x,y,.65,'k')
    hold on
    plot([uData.medianRmse uData.medianRmse],[0 400],'k')
    title(sprintf('PVE Data1/Data2 (%2.2f) | Global (%2.2f) | Voxel-wise (%2.2f)\nRMSE Ratio, Total (%2.2f) | Median (%2.2f)', ...
      uData.totVarExpD, ...
      uData.totVarExpM, ...
      uData.totVarExpMvw, ...
      uData.totrmse,   ...
      uData.medianRmse))
    xlabel('Root mean square error')
    ylabel('Number of occurrences')

  case 'fiberlength'
    % Plot a histogram of the length of the fibers in the connectome.
    % uData = fePlot(fe,'fiber length')
    uData.len       = fefgGet(feGet(fe,'fg acpc'),'length');
    uData.nanmedian = nanmedian(uData.len);
    
    % male a histogram of the length of the fibers
    [y,x] = hist(uData.len,ceil(length(uData.len)/10));
    bar(x,y,2,'r')
    hold on
    plot([uData.nanmedian uData.nanmedian],[0 1400],'k')
    title(sprintf('Total number of fibers: %i',feGet(fe,'n fibers')))
    xlabel('Fiber length in mm')
    ylabel('Number of occurrences')
    
  case 'qualityoffit'
    % Plots the quality of fit of the LiFE model
    % fData = fePlot(fe,'quality of fit');
    uData = fePlotQualityOfFit(fe);
    
  case {'mmatrix','model'}
    % Show an image of the LiFE model/matrix
    % uData = fePlot(fe,'model')
    uData.M = feGet(fe,'M fiber');
    imagesc(uData.M);
    
  case 'dsigmeasured'
    % Plot the series of the measured diffusion signal reordered in such a
    % way that close angle in on the sphere are close together in the
    % series. This call uses the functionality in feGet(fe,dsig measured',
    % varargin)
    % uData = fePlot(fe,'dsig measured')
    % uData = fePlot(fe,'dsig measured',voxelIndex)
    % uData = fePlot(fe,'dsig measured',coords)
    
    % Use routine in dwiPlot to reorder the x-axis of this so that
    % every voxel is plotted
    dsig = feGet(fe,'dsig measured');
    dsig = dsig(1:150);
    
    % From dwiPlot
    bvecs = feGet(fe,'bvecs');
    
    % Prepare the signal in 2D (imag) format:
    bFlat = sphere2flat(bvecs,'xy');
    
    % Interpolating function on the 2D representation
    F = TriScatteredInterp(bFlat,dsig(:));
    
    % Set the (x,y) range and interpolate
    x = linspace(min(bFlat(:,1)),max(bFlat(:,1)),sqrt(feGet(fe,'nbvecs')));
    y = linspace(min(bFlat(:,1)),max(bFlat(:,1)),sqrt(feGet(fe,'nbvecs')));
    [X Y] = meshgrid(x,y);
    est = F(X,Y);
    est(isnan(est)) = 0;  % Make extrapolated 0 rather than NaN
    %est = round(est./2);
    
    inData.x = X;
    inData.y = Y;
    inData.data = est;
    
    uData = dwiSpiralPlot(inData);
    dwiSpiralPlot(inData);

  case 'bvecs'
    % Shows the bvecs on a sphere.
    % uData = dwiPlot(fe,'bvecs')
    uData = dwiPlot(feGet(fe,'dwi'),'bvecs');
    
  case {'connectome','fg','fibers'}    
    % Plot the fiber group in a matlab 3D mesh.
    % uData = fePlot(fe,'fg')
    % We also look for fibers that are 'broken.' And display those as
    % separated fibers.
    uData = feConnectomeDisplay(feSplitLoopFibers( feGet(fe,'fg img')),g);
    %uData = feConnectomeDisplay(feGet(fe,'fg img'),figure);

  otherwise
    error('Unknown plot type %s\n',plotType);
end

% Set Figure properties
if ~notDefined('uData')
  set(g,'userdata',uData)
end

% Set axes properties
set(gca, 'FontSize', 16, ...
  'tickDir','out', ...
  'box','off');

end


%% Scratch
% g = mrvNewGraphWin;
% u = randn(128,128);
% clear uData
% uData.u = u;
% imagesc(u);
% set(g,'userdata',uData);
% 
% foo = get(g,'userdata')

