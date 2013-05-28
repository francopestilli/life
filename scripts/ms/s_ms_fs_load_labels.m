function s_ms_fs_load_labels(subject,hemisphere,recomputelabels)
%
% Load a series of FreeSurfer Labels into an ROI.
%
%
% Written by Franco Pestilli, Stanford University 2013
    
fsSubDir = getenv('SUBJECTS_DIR');

% This is the subject folder. The FreeSurfer folder for the subject. It is
% a folder under $SUBJECTS_DIR
if notDefined('subject')
    subject = 'pestilli_test';
end

% Hemisphere, 'lh' or 'rh'
if notDefined('hemisphere')
    hemisphere = 'lh';
end

% The default label and ROIs dir for this subject.
if notDefined('roiSaveDir')
    roiSaveDir = fullfile(fsSubDir,subject,'label');       
end

% The following are the file name sof the labels we want to load.
% The fullpath to a .label FreeSurfer file.
if notDefined('labels')
    labels = {sprintf('%s.%s.label',hemisphere,'middletemporal'),...
              sprintf('%s.%s.label',hemisphere,'parsopercularis'), ...
              sprintf('%s.%s.label',hemisphere,'parstriangularis'), ...
              sprintf('%s.%s.label',hemisphere,'bankssts')};
          
end
labels = {sprintf('%s.%s.label',hemisphere,'insula'), ...
          sprintf('%s.%s.label',hemisphere,'G_and_S_cingul-Ant')};

% The fullpath to the .nii.gz file that will be saved out. With NO .nii.gz 
% extension.
if notDefined('niftiRoiName')
    niftiRoiDir = fullfile(roiSaveDir,'nifti');
    if ~isdir(niftiRoiDir)
        mkdir(niftiRoiDir)
    end
    for ir = 1:length(labels)
        niftiRoiName{ir} = fullfile(niftiRoiDir,labels{ir});
    end  
end

if notDefined('recomputelabels')
    recomputelabels = 1;
end

if recomputelabels
    % % In case not all the labels were create dtis is an example of how to
    % % create labesl from a segementation
    annotation     = 'aparc';
    annotationFile = fullfile(fsSubDir,subject,'label',sprintf('%s.%s.annot',hemisphere,annotation));
    fs_annotationToLabelFiles(subject,annotationFile,hemisphere)
end

for ir = 1:length(labels)
    fs_labelFileToNiftiRoi(subject,labels{ir},niftiRoiName{ir},hemisphere);
end

end % end main function
