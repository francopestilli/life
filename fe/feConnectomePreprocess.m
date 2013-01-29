function [fg fibersToKeep] = feConnectomePreprocess(fg, volumeRoi, cortex, cortVals, minCorticalDist, maxVolDist, min_length)
% Preprocess a connectome. Always remove fibers that are too short and
% optionally remove fibers that:
%    1. Do not terminate in the gray matter.
%    2. Exit a predefiend volume of interest (The connectome ROI).
%
% Please note that this function is meant to be used to preprocess a
% whole-brain tracktogaphy connectome. Using it for smaller connectomes can
% introduce major errors in the fiber selection process.
%
% [fgOut fibersToKeep] = feConnectomePreprocess(fg, cortex, minCorticalDist, ...
%                        minCorticalDist, volumeRoi, maxVolDist, min_length)
%
% INPTUS:
%    fg              - The connectome to preprocess. THis is generally a
%                      larger fiber group as defined in mrDiffusion.
%    cortex          - This is a segmentation to use to identify the cortex. 
%                      It can be passed in either as:
%                      * A path to a NIFTI classification file, like the one 
%                        returned by freesurfer.
%                      * A mrDiffusion ROI.
%                      * A NIFTI file.
%    cortVals        - A vector of values that define cortex within the nifti image.
%                      Default are the values in a freesurfer class file (5 6)
%    minCorticalDist - Distance criteria for what counts as terminating in the
%                      cortex. The default is 2mm meaning that fibers must terminate
%                      within 2mm of a cortex voxel will be retained.
%    volumeRoi       - An cortexRoi in the same coordinate space as the fiber group.  
%                      This ROI should be the volumeRoi containing the final fiber
%                      group. The fibers will be clipped upon exiting this volumeRoi.
%                      All further LiFE analyses will be confined to this region.
%    maxVolDist      - A fiber node will be considered outside of the cortexRoi when it
%                      is more than maxVolDist mm away from any of the coordinates
%                      in the volumeRoi cortexRoi.
%    min_length      - Minimum fiber length. Fibers with fewer nodes (number of
%                      xyz coordinates) will be discarded.
%
% OUTPUTS:
%    fgOut           - The preprocessed fiber group.
%    fibersToKeep    - A vector of one's and zero's defining which fibers from the 
%                      original group were kept after preprocessing.
%
% EXAMPLE:
%    cortex = '/biac4/wandell/biac2/wandell2/data/WMDevo/adult/anatomy/107_JW/t1_class.nii.gz'
%    fg     = mtrImportFibers('/biac4/wandell/biac2/wandell2/data/WMDevo/adult/107_JW/DTI/20111027_1340/dti96ls/fibers/L_SLFFG.pdb')
% 
%    % We generate a spehrical cortexRoi around a portion of the fiber group.
%    volumeRoi        = dtiNewRoi;
%    volumeRoi.coords = dtiBuildSphereCoords(mean(fg.fibers{1}'),20);
%    cortVals = 1;
%    fgOut                = feConnectomePreprocess(fg, volumeRoi, cortex, cortVals)
%
% Jason & Franco (c) 2012 Stanford Vista Team.

%% Check arguments
if notDefined('fg'), error('[%s] Please supply a fiber group', mfilename);end

% Fiber length is defined in nodes (e.g., number of xyz coordinates in a
% fiber).
if notDefined('min_length'),  min_length = 20; end

% Minimum distance of a node from the cortex. Default to 1mm.
if notDefined('maxVolDist'),  maxVolDist = 1;end

initial_f_num = numel(fg.fibers);

%% Find long fibers
tic, fprintf('[%s] Finding long fibers...',mfilename);
long_fibers = feFindLongFibers(fg.fibers,min_length);
fprintf('process completed in %2.3fs.\n',toc);

%% Find fibers with both fiber endpoints within minCorticalDist from the cortex
if ~notDefined('cortex')
    tic, fprintf('[%s] Finding fibers intersecting the cortex...',mfilename);
    % Values identifying the left and rigth cortex in a classification file.
    % Default to freesurfer values.
    if notDefined('cortVals'), cortVals = [5 6];end
    % Default is to discard fibers that are more than 2mm from cortex.
    if notDefined('minCorticalDist'), minCorticalDist = 3;end
    cortexRoi               = feBuildCorticalRoi(cortex, cortVals);
    [~, ~, cortical_fibers] = dtiIntersectFibersWithRoi([], ...
                              {'and', 'both_endpoints'},    ...
                              minCorticalDist, cortexRoi, fg);
    fprintf('process completed in %2.3fs.\n',toc);
else
  cortical_fibers = true(length(fg.fibers)); 
end

%% Extract the fibers that pass the requested criteria.
tic, fprintf('[%s] Extracting long fibers that intersect the cortex...',mfilename);
fibersToKeep = (cortical_fibers & long_fibers);
fg           = fgExtract(fg,find(fibersToKeep), 'keep');
fprintf('process compelted in %2.3fs.\n', toc);

%% Clip fibers to be constrained within the volumeRoi.
tic, fprintf('[%s] Clipping fibers that leave the volume ROI...\n',mfilename);
fg.fibers = feClipFibersToVolume(fg.fibers,volumeRoi.coords,maxVolDist);
fprintf('process completed in %2.3fhours\n',toc/60/60);

fprintf('[%s] Results:\n- %i fibers deleted out of %i initial fibers (%i%%).\n- %i short fibers,\n- %i non ending in cortex.\n', ...
    mfilename,length(find(~fibersToKeep)),initial_f_num,round(100*(length(find(~fibersToKeep))/initial_f_num)),...
    length(find(~long_fibers)),length(find(~cortical_fibers)));

end % End main function.

