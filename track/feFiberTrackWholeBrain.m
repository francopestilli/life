function [fg outName] = feFiberTrackWholeBrain(dt6,trackOpts,outName,hem,save,fileFormat,wmSeedRegion)
% 
% function [fg outName] = mctFiberTrackWholeBrain([dt6],[trackOpts=dtiFiberTrackOpts],...
%               outName=['FG_ ' hem],[hem='both'],[save=1]),[fileFormat='pdb'],[wmSeedRegion])
% 
% This function will return a whole-brain group of fibers given some dt6
% file and some tracking parameters. User can specify which hemisphere to
% return with the 'hem' input argument. By default the whole-brain is
% tracked using a moderately thresholded FA mask of 0.35 as an ROI. Most
% default tracking parameters from dtiFiberTrackOpts are used if the user
% does not pass them in (see example usage). To change the number of fibers
% returned the user should edit the trackOpts by following the example
% below. Note that the default trackOpts return ~200K (depending on the
% data) and takes ~4-5min to run, thus, in this function - if the user does
% not pass in trackOpts - seedVoxelOffsets is set to 0, to return a
% reasonable number of pathways (~40k). If the user does pass in trackOpts
% we assume they know what they are asking for. 
% 
% INPUTS:
%       dt6        - File obtained from dtiInit
%       trackOpts  - a structure containing the tracking params. Set with
%                    dtiFiberTrackOpts. The default values are:
%                        opts.stepSizeMm = 1; % mm
%                        opts.faThresh = 0.20;
%                        opts.lengthThreshMm = 20;
%                        opts.angleThresh = 30;
%                        opts.wPuncture = 0.2;
%                        opts.whichAlgorithm = 1;
%                        opts.whichInterp = 1;
%                        opts.seedVoxelOffsets = [-0.25  0.25];
%                        opts.offsetJitter = 0;
%       outName    - Name of the output fibers as they will be saved.
%       hem        - The hemisphere for which you want paths tracked.
%                    Default is 'both'. 'left' or 'right' optional. 
%       fileFormat - 'pdb' or 'mat' to save the output fibers.
%       save       - 1 = save fibers, 0 = don't save the fibers (useful if
%                    you just want to return the struct).
%       wmSeedRegion     - A file of white matter coordinates. It will be used as ROI 
%                    for tracking isntead of creating an ROI out of a FA-thresholded 
%                    mask. The mask is assumed to ba at the resolution of
%                    the dwi data and containing only 1's and 0's. Where 1
%                    indicates white matter. (see feMakeWhiteMatterMask.m
%                    for an example of how to create a WM mask.)
%
% OUTPUTS:
%       fg         - structure containing the fibers in fg.fibers.
% 
% WEB RESOURCES:
%       mrvBrowseSVN('mctFiberTrackWholeBrain');
% 
% EXAMPLE:
%       trackOpts = dtiFiberTrackOpts('faThresh',.30,'seedVoxelOffsets',[0.25]);
%       dt6       = 'dt6.mat';
%       outName   = 'WholeBrainFG';
%       hem       = 'both';
%       fg = dtiFibersTrackWholeBrain(dt6,trackOpts,outName,hem,fileFormat,save)
% 
% NOTES: Based on dtiFiberTrackWholeBrain
%
% Franco (C) Stanford VISTASOFT Team 2012 

% Check inputs
fg = [];
if notDefined('dt6')
    dt6 = mrvSelectFile('r','*.mat','Please select a dt6.mat file.');
end

% If the user does not pass in trackOpts we set them to default values
% (except for seedVoxelOffsets, which is set to 0 to return a
% reasonable number of pathways (~40k)). 
if notDefined('trackOpts');  trackOpts  = dtiFiberTrackOpts('seedVoxelOffsets', 0) ; end
if notDefined('hem');        hem        = 'both'; end
if notDefined('fileFormat'); fileFormat = 'pdb'; end
if notDefined('outName');    outName    = ['fe_whole_brain_' hem 'Hem']; end
if notDefined('save');       save       = 1; end

% We don't want the full path in the fiber.name field. We will rebuild
% later in the case that the user want's the fibers saved.
[p f] = fileparts(outName); outName = f; 

%% Switch on possible inputs for hem and fileFormat
switch lower(hem)
    case {'left' ,'l'}; hem = 1;
    case {'right','r'}; hem = 2;
    case {'both' ,'b'}; hem = 0;
end

switch lower(fileFormat)
    case {'pdb','.pdb','p'}; fileFormat = 1;
    case {'mat','.mat','m'}; fileFormat = 0;
end            

% Load dt6 
if     ischar(dt6), dt = dtiLoadDt6(dt6);
elseif iscell(dt6), dt = dt6;
else                keyboard;
end
  
% If a white-matter mask was not passed in we guess a reasonable
% white-matter mas using FA thresholding. This is not quite correct because
% reagions of corssing fibers might not pass the criteria and enter as seed
% for tracking. Btu this guess has worked reasonably well in the past,
% because the mask here is only a seed mask and the actual FA thresholding
% is done inside the c-code of dtiFiberTracker.cxx
if notDefined('wmSeedRegion')
    % Compute FA and create a brain mask ROI
    %
    % Compute the FA and normalize
    fa = dtiComputeFA(dt.dt6);
    fa( fa > 1 ) = 1;
    fa( fa < 0 ) = 0;
    
    % Using a reasonable FA as a threshold - if you use opts.faThresh (.15) for
    % this purpose you get ~1 million paths back, which takes way too long to
    % track. Instead we only use voxels with an FA of at least .35 as seed
    % points. This gives us ~200-400K paths with default trackOpts.
    faThresh = 0.35;
    mask 	 = fa >= faThresh;
else
    % The white matter mask is assumed to have only 0's and 1's inside.
    % Where 1's indicate white-matter
   ni = niftiRead( wmSeedRegion );
   mask = ni.data; 
   clear ni
end

% Find the voxels that we will use as sed for tracking
[x,y,z]  = ind2sub(size(mask), find(mask));

% Create a new whole-brain ROI from the mask
roiAll        = dtiNewRoi('all');
roiAll.coords = mrAnatXformCoords(dt.xformToAcpc, [x,y,z]);

if ~notDefined('wmSeedRegion')
    % Smooth the ROI and fill holes
    roiAll = dtiRoiClean(roiAll,3,{'fillHoles'});
end

% Generate ROI coords for each case
% 
% Whole brain (hem = 'both')
if hem == 0; disp('Tracking Whole-Brain Fibers ...');
    coords = roiAll.coords; 
end

% Left hemisphere (hem = 'left')
if hem == 1; disp('Tracking Left Hemisphere Fibers ...');
    roiLeft = dtiRoiClip(roiAll, [1 80]);
    coords  = roiLeft.coords;
end

% Left hemisphere (hem = 'right')
if hem == 2; disp('Tracking Right Hemisphere Fibers ...');
    roiRight = dtiRoiClip(roiAll, [-80 -1]);
    coords   = roiRight.coords;
end

if ~(exist(fullfile(p,[outName,'.mat']),'file') == 2)
  % Track and clean the fibers
  % fg = dtiFiberTrack2(dt.dt6, coords, 1, 0, dt.mmPerVoxel, dt.xformToAcpc, outName, trackOpts);
  
  % Track and clean the fibers - using STT
  fg = dtiFiberTrack(dt.dt6, coords, dt.mmPerVoxel, dt.xformToAcpc, outName, trackOpts);

  % Clean the fibers. This will ensure there is only one point that crosses
  % the midline.
  % fg = dtiCleanFibers(fg); % I commented this out because the cleaning
  % resamples the fibers at 1mm.
  
  % Save out the fibers using the selected fileFormat (mat or pdb)
  if save == 1
    % Rebuild outName in the case that the user passed in a full path
    % originally.
    outName = fullfile(p,outName);
    if(fileFormat == 0); dtiWriteFiberGroup(fg, outName); end
    if(fileFormat == 1); mtrExportFibers(fg, outName);    end
    fprintf('Saved %s. \n', outName);
  end
else
  outName = fullfile(p,outName);
  fprintf('\nFound tracking file %s. Not tracking.\n', outName);
end
  
return
