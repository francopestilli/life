function s_fe_import_roi_from_mrvista
%
% Example script for importing ROI from mrVISTA to mrDiffusion.
%
%
% Franco (c) 2012 Stanford VISTA Team

% Choose a data set
dataRootPath = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
subfolders   = fullfile('150dirs_b1000_1');
baseDir      = fullfile(dataRootPath,subfolders);
dtFile       = fullfile(baseDir,'dt6.mat');
t1File       = fullfile(dataRootPath,'t1','t1.nii.gz');
saveDir      = fullfile(baseDir,'ROIs');
% ROIs: /biac3/wandell7/data/Retinotopy/Pestilli/mrvSession1/Gray/ROIs/

% Now import some ROIs by selecting them from the GUI.
dtiRoiXformMrVistaVolRoi(dtFile,[],t1File,saveDir,'nifti');
