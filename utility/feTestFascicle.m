function [feWithoutFas, feWithFas, commonCoords] = feTestFascicle(fe,fascicleIndices,displayFibers)
%
% [feWithoutFas, feWithFas, commonCoords] = feTestFascicle(fe,fascicleIndices,[displayFibers])
%
% This file tests the hypothesis of the importance of a fascicle in a
% volume of white matter.
%
% Franco (c) Stanford Vista Team 2013

if notDefined('displayFibers'), displayFibers = 0;end

% Remove the fibers of the fascicle from the fe.
feWithoutFas = feConnectomeReduceFibers(fe, ~fascicleIndices );

% Extract the fascicle out of the fiber group.
fas = fgExtract(feGet(fe,'fibers img'), fascicleIndices, 'keep' );

% Get the cordinates of the fascicle we just deleted. These coordinates are
% contained in the full connectome. We want to fin the indices in the M
% matrix to the the voxels of the fascicle and keep only those voxels.
fasCoords    = fefgGet(fas,'unique image coords');
allCoords    = feGet(fe,'roi coords');
commonCoords = ismember(allCoords, fasCoords,  'rows');

% Now: commonCoords contains the indices of the voxels to keep from feFas,
% these voxels are part of the connectome fe. So now we want to delete
% all the rest of the voxels from feWithoutFasFP and keep only the voxels where
% the fascicle in feFP goes through. At the same time we keep all the
% fibers in the connectome feWithoutFasFP and test the quality of fit of the
% connectome at the location where the feFN fascicle went throguh.
feWithoutFas = feConnectomeReduceVoxels(feWithoutFas,find(commonCoords));
feWithFas    = feConnectomeReduceVoxels(fe,          find(commonCoords));

% Refit and install the fit.
fprintf('\n[%s] Fitting connectome WITHOUT the fascicle.',mfilename)
feWithoutFas = feSet(feWithoutFas,'fit', feFitModel(feGet(feWithoutFas,'Mfiber'),feGet(feWithoutFas,'dsigdemeaned'),'sgdnn'));
fprintf('\n[%s] Fitting connectome WITH the fascicle.',mfilename)
feWithFas    = feSet(feWithFas,   'fit', feFitModel(feGet(feWithFas,'Mfiber'),   feGet(feWithFas,'dsigdemeaned'),   'sgdnn'));

% The following is test code to show where the coordinates of the facicle
% that were removed land inside the connectoem. Also I show in gray hte
% connectome WITHOUT the fascicle and in red the connectome WITH the
% fascicle. Where there is only red in the connectoem that is where the
% fascicle was removed but we are attemtping to explain the variance in the
% data.
if displayFibers
    % Show the coordinates to see if we are in the right spot
    mrvNewGraphWin('Coordinate check');
    plot3(allCoords(commonCoords,1),allCoords(commonCoords,2),allCoords(commonCoords,3),'ko','MarkerFaceColor','k','MarkerSize',8);
    hold on;
    plot3(fasCoords(:,1),fasCoords(:,2),fasCoords(:,3),'ro','MarkerFaceColor','r','MarkerSize',3);
    axis equal
    %view(3,79)
    view(-23,-23);
    
    % Now display the fascicle removed and the connectome without the fascicle
    % in one figure with different colors.
    %
    feConnectomeDisplay(feSplitLoopFibers(feGet(feWithoutFas,'fibers img')),figure);
    hold on
    %feConnectomeDisplay(feSplitLoopFibers(feGet(feWithFas,'fibers img')),gcf, [.95 .1 .1])   
    feConnectomeDisplay(feSplitLoopFibers(fas),gcf,[.95 .1 .1])
    view(-23,-23);
    h= camlight;

    feConnectomeDisplay(feSplitLoopFibers(feGet(feWithoutFas,'fibers img')),figure);
    hold on
    plot3(allCoords(commonCoords,1),allCoords(commonCoords,2),allCoords(commonCoords,3),'go','MarkerFaceColor','g','MarkerSize',12);
    view(-23,-23);
    h= camlight;       
end
 
end