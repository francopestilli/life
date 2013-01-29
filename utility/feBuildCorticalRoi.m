function cortexRoi = feBuildCorticalRoi(cortex, cortVals)
%
% Generate an ROI from a NIFTI classification file (e.g., a freesurfer-class).
% 
%    cortexRoi = feBuildCorticalRoi(cortex, cortVals)
% 
% INPUTS: 
%        cortex  - A volume to use to identify the cortex.
%                  It can be passed in either as a:
%                  * Path to a NIFTI classification file, like the one returned by
%                    freesurfer.
%                  * mrDiffusion ROI
%                  * NIFTI file.
%        cortVals - A vector of values identifying the cortical surface in
%                   the classification file. Default to freesurfer values
%                   (5, 6).
%
% OUTPUTS:
%         cortexRoi - A mrDiffusion ROI defining the cortical surface in
%         x,y,z coordinates.
%
% SEE ALSO: feConnectomePreprocess.m
%
% Franco (c) 2012 Stanford Vista Team.

%% Check inputs
if ~exist('cortex','var') || isempty(cortex)  % No cortex was passed, return an error.
  fprintf('[%s] Please supply a ''cortex'' volumeRoi as a NIFTI file.',mfilename)
elseif ischar(cortex) && exist(cortex,'file') % A path to a nifti was passed.
  % Do nothing.
elseif isstruct(cortex) && isfield(cortex,'data') % Nifti image was passed in
  % Down we use dtiRoiFromNifti to create a cortexRoi. The function only handles
  % paths to nifti files. So here we take the path from the nifti file name
  % and reload it, later. Not optimized for speed but reliable.
  cortex = cortex.fname; 
elseif istruct(cortex) && isfield(cortex,'coords')
  % If an cortexRoi containing a list of cortex coordinates was passed in
  % (rather than an image) use that
  cortex = [];  cortexRoi = cortex; clear cortex
else error('[%s] ''cortex'' must be either NIFTI or a path to a NIFTI.', mfilename)
end

% Values identifying the left and rigth cortex in a classification file.
% Default to freesurfer values.
if ~exist('cortVals','var') || isempty(cortVals), cortVals = [5 6];end

%% Generate a cortical ROI from the NIFTI file on disk
if ~isempty(cortex)
  cortexRoi = dtiRoiFromNifti(cortex, cortVals,'cortex_mask','mat',[],false);
else
  error('[%s] Could not geenrate a cortical ROI, please check inputs.',mfilename)
end

end % End main function
