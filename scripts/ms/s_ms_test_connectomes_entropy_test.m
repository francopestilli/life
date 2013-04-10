function s_ms_test_connectomes_entropy_test(trackingType,lmax,diffusionModelParams)
%
%Preprocess a connectome (fg) to be constrained within a
% region of interest and within the cortex.
%
% Fit life to it. Returns a series of maps and saves out the results in the
% folders for the ms on azure
%
% For this example we use the volume defined by the connectome to constrain
% the fibers.
%
% See also:
%
% Written by Franco Pestilli (c) 2012 Stanford University, VISTA team.

if notDefined('trackingType'),trackingType = 'probabilistic';end
if notDefined('lmax'),        lmax=[6];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('dataType'), dataType='150dirs';end

%thisRoi = {'fp_right_sphere_blue_UP_DOWN','fp_right_sphere_green_FRONT_BACK','fp_right_sphere_red_LEFT_RIGHT'};
thisRoi = {'fp_corpus_collosum'};

% Overwrite files
clobber = 1;

%% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

for iroi = 1:length(thisRoi)
%% Start loading the data
% Basic directories
switch dataType
    case {'150dirs','2mm'}
        %% DWI data
        dataRootPath  = fullfile('/biac2/wandell2/data/diffusion/pestilli/20110922_1125');
        subfolders    = fullfile('150dirs_b2000_1');
        baseDir       = fullfile(dataRootPath,subfolders);
        dtFile        = fullfile(baseDir,'dt6.mat');
        dwiFile       = fullfile(dataRootPath,'raw','0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin.nii.gz');
        dwiFileRepeat = fullfile(dataRootPath,'raw','0007_01_DTI_2mm_150dir_2x_b2000_aligned_trilin.nii.gz');
        t1File        = fullfile(dataRootPath,'t1','t1.nii.gz');
        
    
    %% ROIs connectomes and saved paths
    projectDir           = '/azure/scr1/frk/150dirs_b1000_b2000_b4000';
    connectSubfolders    = {'life_mrtrix_rep1'};
    testVol              = sprintf('/azure/scr1/frk/150dirs_b1000_b2000_b4000/ROIs_restric_connectomes/%s.mat',thisRoi{iroi});
    saveDir              = fullfile(projectDir,sprintf('results_%s',thisRoi{iroi}));
    mkdir(saveDir);
    
    % Get the connectome names
    switch trackingType
       case {'probabilistic','p'}
        connectomeFile = { ...
          sprintf( '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax%i_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb',lmax),...
          };
       otherwise
        keyboard
    end
  
  otherwise
    keyboard
end

% Max distance in mm from the ROI edges.
maxVolDist = 2;

% Load the ROI defining the volume we will use to test the connectome.
volRoi = dtiReadRoi(testVol);
for irep = 1:3
  for ii = 1:length(connectomeFile)
    fprintf('Processing: \n %s \n ======================================== \n\n',connectomeFile{ii})
    
    %% Set up file names
    % Find an identificative name for the connectome that is short
    % enough:
    cName = [connectomeFile{ii}(1:57),'_',connectomeFile{ii}(end-16:end-4)];
    feSaveDir         = fullfile(saveDir,connectSubfolders{irep},'fe_structures');
    feSaveName        = sprintf('%s_diffModAx%sRd%s_%s',cName,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)));
    switch dataType
      case {'150dirs','2mm'}
        
        % Select the dwi file
        if ii == 1
          dataRootPath  = fullfile('/biac2/wandell2/data/diffusion/pestilli/20110922_1125');
          subfolders    = fullfile('150dirs_b2000_1');
          baseDir       = fullfile(dataRootPath,subfolders);
          dtFile        = fullfile(baseDir,'dt6.mat');
          dwiFile       = fullfile(dataRootPath,'raw','0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin.nii.gz');
          dwiFileRepeat = fullfile(dataRootPath,'raw','0007_01_DTI_2mm_150dir_2x_b2000_aligned_trilin.nii.gz');
               else
          keyboard
        end
    end
    
    if ~(exist(fullfile(feSaveDir,[feSaveName,feSaveName,'.mat']),'file') == 2) || clobber
      fgSaveName        = fullfile(saveDir,connectSubfolders{irep},'fg',feSaveName);
  
      % Create the necessary folders if they were not there already
      foldersToCheck = {fgSaveName,fullfile(feSaveDir,'tem.mat')};
      checkFolders(foldersToCheck);
      
      %% Load a full-brain Connectome that was preprocessed with...
      wholeBrainPath = fullfile(projectDir,connectSubfolders{irep},connectomeFile{ii});
      fprintf('Loading the whole brain connectome...\n%s\n',wholeBrainPath)
      fg = fgRead(wholeBrainPath);
      
      %% Clip the connectome to be constrained within the volumeRoi.
      tic, fprintf('Clipping Connectome fibers that leave the volume ROI...\n');
      fg = feClipFibersToVolume(fg,volRoi.coords,maxVolDist);
      fprintf('process completed in %2.3fhours\n',toc/60/60);
      
      %% Build a model of the connectome.
      fe = feConnectomeInit(dwiFile,dtFile,fg,feSaveName,feSaveDir,dwiFileRepeat,t1File,diffusionModelParams);
      feConnectomeSave(fe,feSaveName);
      fgWrite(fg, fgSaveName,'pdb');
      clear fg
      
      %% Fit the model with global weights.
      fefit = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsig demeaned'),'sgdnn');
      fe    = feSet(fe,'fit',fefit);
      feConnectomeSave(fe,feSaveName);
      
          clear fe
      fprintf('DONE processing: \n%s\n ======================================== \n\n',connectomeFile{ii})
    else
      fprintf('FOUND FILE NOT PROCESSING: \n%s\n ======================================== \n\n\n',connectomeFile{ii})
    end
  end
end
end
% Handling paralle processing.
if ~poolwasopen, matlabpool close; end

end % Main function

%-----------------------%
function fold = checkFolders(foldersToCheck)
% Make sure the folders exist otherwise create them:
for ff = 1:length(foldersToCheck)
  [fold,fil,ext] = fileparts(foldersToCheck{ff});
  if ~( exist(fold,'dir') == 7)
    mkdir(fold)
  end
end
end