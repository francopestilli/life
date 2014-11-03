% life_test
feOpenLocalCluster;

%% Build the file names for the diffusion data, the anatomical MRI.
dwiFile       = fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz');
dwiFileRepeat = fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan2_subject1_b2000_150dirs_stanford.nii.gz');
t1File        = fullfile(lifeDemoDataPath('anatomy'),  'life_demo_anatomy_t1w_stanford.nii.gz');

%% 
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
                'life_demo_mrtrix_csd_lmax10_probabilistic.mat');

fg = fgRead(fgFileName);
small_fg = fgCreate('name', ['small_' fg.name], 'fibers', fg.fibers(1:2));
small_fg.pathwayInfo = fg.pathwayInfo(1:2);
fgWrite(small_fg, fullfile(lifeDemoDataPath('diffusion'), 'small_fg.mat'));

%% 
feFileName    = 'life_build_model_demo_CSD_PROB_small';

%% 
fe = feConnectomeInit(dwiFile, small_fg, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);
%fe = feConnectomeInit(dwiFile, fgFileName, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);

%%  
fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

%% 
rmse   = feGet(fe,'vox rmse');
