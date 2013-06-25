function feSaveBoldRoiToNifti(mrsDir, ROIs)
%
% Save a set of mrVista ROIs as a nifti-1 files at the resolution of the T1 anatomy file.
% 
% See also roiSaveForItkGray in vistasoft
%
% Written by Franco Pestilli (c) Vistasoft, Stanford University 2013
mrGlobals;

% Select the mrVista session
if notDefined('mrsDir')
   mrsDir = '/biac4/wandell/biac2/wandell6/data/frk/prf/fp20110912/mrvSession1/';
end

% ROI names
if notDefined('ROIs'), ROIs = {'RTO1-HT','RTO2-HT','RIPS0','RIPS1','RIPS2','RIPS3','LTO1','LTO2-HT','LIPS0','LIPS1','LIPS2','LIPS3'}; end

% Set up a directory for saving the NIFTI ROIs
if notDefined('roisDir')
   roisDir = msPaths('mtrois');
   mkdir(roisDir)
end

% Voxels contained in the ROI are marked as a 1
if notDefined('roiIndex'), roiIndex = 1; end

% Move to the folder containing the mrSession 
cd(fullfile(mrsDir))

% Start hadling a mrVista view:
vw = initHiddenGray; % Open a hidden Gray view for mrVista
vw = loadAnat(vw);   % Load the anatomy
vw = loadROI(vw, ROIs,[], [], [], true);% Load the ROIs

% Create nifti's out of each ROI
for iRoi     = 1:length(ROIs)
    % Build an ROI file
    roiFileName = fullfile(roisDir, [ROIs{iRoi} '.nii.gz']);
    
    % Set the current ROI in the view
    vw = viewSet(vw, 'Current ROI', iRoi);
    
    % Save out the current ROI as a nifti aligned to the T1 Anatomy:
    roiSaveForItkGray(vw, roiFileName,roiIndex)
end

end % End main function

