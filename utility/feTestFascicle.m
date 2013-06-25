function [feWithoutFas,  feWithFas, connectivity, newFascicleIndices, indicesFibersKept, commonCoords] = ...
          feTestFascicle(feWithFas,                  fascicleIndices, refitConnectome, displayFibers)
%
% [feWithoutFas, feWithFas, connectivity, newFascicleIndices, indicesFibersKept, commonCoords] = ...
%  feTestFascicle(feWithFas,fascicleIndices,[refitConnectome],[displayFibers])
%
% This file tests the hypothesis of the importance of a fascicle in a
% volume of white matter.
%
% Written by Franco Pestilli (c) Vista Team, Stanford University 2013

if notDefined('displayFibers'),   displayFibers = 0;end
if notDefined('refitConnectome'), refitConnectome=1;end

% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

% Remove the fibers of the fascicle from the feWithFas structure.
feWithoutFas = feConnectomeReduceFibers(feWithFas, ~fascicleIndices );

% Extract the fascicle out of the fiber group.
fas = fgExtract(feGet(feWithFas,'fibers img'), fascicleIndices, 'keep' );

% Get the cordinates of the fascicle we just deleted. These coordinates are
% contained in the full connectome. We want to fin the indices in the M
% matrix to the the voxels of the fascicle and keep only those voxels.
fasCoords    = fefgGet(fas,'unique image coords');
allCoords    = feGet(feWithFas,   'roi coords');
commonCoords = ismember(allCoords, fasCoords,  'rows');

% Now: commonCoords contains the indices of the voxels to keep from feFas,
% these voxels are part of the connectome feWithFas. So now we want to delete all
% the rest of the voxels from feWithoutFasFP and keep only the voxels where
% the fascicle in feFP goes through.
% 
% At the same time we keep all the fibers in the connectome that still pass
% throught the left voxels in the ROI. We will use this subset of voxels
% and fibers to test the quality of fit of the connectome at the location
% where the feFN fascicle went through.
feWithoutFas = feConnectomeReduceVoxels(feWithoutFas,find(commonCoords));
[feWithFas, indicesFibersKept]    = feConnectomeReduceVoxels(feWithFas, find(commonCoords));

% Here we return the indices of the fascicle in the newly resized
% connectome.
newFascicleIndices = ismember(find(indicesFibersKept),...
                              find(fascicleIndices),'rows');

% Fit the connectome again with the left fibers.
if refitConnectome
    % Refit and install the fit.
    fprintf('\n[%s] Fitting connectome WITHOUT the fascicle.',mfilename)
    feWithoutFas = feSet(feWithoutFas,'fit', ...
        feFitModel(feGet(feWithoutFas,'Mfiber'),feGet(feWithoutFas,'dsigdemeaned'), ...
        'sgdnn'));
    
    fprintf('\n[%s] Fitting connectome WITH the fascicle.',mfilename)
    feWithFas    = feSet(feWithFas,   'fit', ...
        feFitModel(feGet(feWithFas,'Mfiber'),   feGet(feWithFas,'dsigdemeaned'), ...
        'sgdnn'));
end

% Now compute the connectivity measure of the voxel.
% This is the sum of the weights of the fascile divided by the sum of the
% weights of all the fibers in the same volume of white-matter.
%
% Find the weights for all the fibers that go through the voxels of
% this connection.
connectivity.wall = feGet(feWithFas,'fiber weights');
connectivity.wfas = connectivity.wall(newFascicleIndices);
connectivity.wnfas = connectivity.wall(~newFascicleIndices);

% Compute some measures of strength of the connection represented by the
% fascicle, by comparing the fascicle weights with the weights of all the
% rest of the fascicles going through the same voxels.
connectivity.strength(1) =    sum(connectivity.wfas) /  sum(connectivity.wnfas);
connectivity.strength(2) =   mean(log10(connectivity.wfas)) / mean(log10(connectivity.wnfas));
connectivity.strength(3) = median(log10(connectivity.wfas))/median(log10(connectivity.wnfas));

% The following is test code to show where the coordinates of the facicle
% that were removed land inside the connectoem. Also I show in gray hte
% connectome WITHOUT the fascicle and in red the connectome WITH the
% fascicle. Where there is only red in the connectoem that is where the
% fascicle was removed but we are attemtping to explain the variance in the
% data.
if displayFibers
    % Show the coordinates to see if we are in the right spot
    mrvNewGraphWin('Coordinate check');
    plot3(allCoords(commonCoords,1),allCoords(commonCoords,2), ...
          allCoords(commonCoords,3),'ko','MarkerFaceColor','k','MarkerSize',8);
    hold on;
    plot3(fasCoords(:,1),fasCoords(:,2),fasCoords(:,3),'ro', ...
          'MarkerFaceColor','r','MarkerSize',3);
    axis equal
    %view(3,79)
    view(-23,-23);
    
    % Now display the fascicle removed and the connectome without the fascicle
    % in one figure with different colors.
    %
    feConnectomeDisplay(feSplitLoopFibers(feGet(feWithoutFas,'fibers img')),figure);
    hold on
    %feConnectomeDisplay(feSplitLoopFibers(feGet(feWithFas,'fibers img')),gcf, [.95 .1 .1])   
    feConnectomeDisplay(feSplitLoopFibers(fas),gcf,[.95 .1 .1]);
    view(-23,-23);
    h= camlight;

    feConnectomeDisplay(feSplitLoopFibers(feGet(feWithoutFas, ...
                        'fibers img')),figure);
    hold on
    plot3(allCoords(commonCoords,1),allCoords(commonCoords,2), ...
        allCoords(commonCoords,3),'go','MarkerFaceColor','g','MarkerSize',12);
    view(-23,-23);
    h= camlight;       
end

if ~poolwasopen, matlabpool close; end

end