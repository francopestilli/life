function s_ms_cortical_connection()
%
% Extrat Broca's area and the middle temporal-gyrus area from a FreeSurfer segementation
%
%
%
%
% Written by Franco Pestilli (c) 
% Stanford University 2013 Vistasoft
fsDir        = getenv('SUBJECTS_DIR');
subject      = 'pestilli_test';
hemisphere    = 'lh'; 
annotation   = 'aparc.a2009s'; % Or aparc.annot (These segmentations are from different papers, see FS online)
segmentation = 'aparc.a2009s+aseg.mgz';

% Create a nifti ROI from a segmenttion file.
%
% Load the brain areas segmented in FreeSurfer
annotationFile = fullfile(fsDir,subject,'label',annotation);
fs_annotationToLabelFiles(subject,annotation)

% Create nifti ROi from the label file that we want:
%
% These ROIs contain Wernicke's area:
superior = {'G_temp_sup-Plan_tempo','S_temporal_transverse',...
            'G_temp_sup-Lateral','S_temporal_sup'};
middle   = {'G_temporal_middle','S_temporal_inf', ...
            'G_temporal_inf'};

% Now we can create Nifti ROIs from all the labels
% Combine them in mrDiffusion or FreeSurfer?
this_roi = 'G_temporal_inf';
labelFileName = sprintf('%s.%s.label',hemisphere,this_roi);
niftiRoiName  = sprintf('%s_%s_label',hemisphere,this_roi);
[niftiRoiName, niftiRoi] = fs_labelFileToNiftiRoi(subject,labelFileName,niftiRoiName);
                       
% Create a nifti ROI from the BA44 and BA45 (Broca's area) provided by
% Amunts
broca = {'BA44','BA45'};
hemisphere = 'lh';
labelFileName = sprintf('%s.%s.label',hemisphere,broca{1});
niftiRoiName  = sprintf('%s_%s_label',hemisphere,broca{2});
[niftiRoiName, niftiRoi] = fs_labelFileToNiftiRoi(subject,labelFileName, niftiRoiName);

