function fe = s_fe_fit_whole_brain
%
% Illustrate how to open up data and run a small linear fascicle evaluation
% (LIFE).
% 
% s_fe_fit_whole_brain
%
% We should start making fePlot(), a gateway routine that takes an fe
% structure and some other arguments and generates visualizations of
% interesting quantities.
%
%  [uData, g] = fePlot(fe,plotType,varargin);
%
%
% Franco (C) 2012 Stanford VISTA team.

dataRootPath = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
subfolders   = fullfile('150dirs_b1000_1');
baseDir      = fullfile(dataRootPath,subfolders);
dtFile       = fullfile(baseDir,'dt6.mat');
dwiFile      = fullfile(dataRootPath,'raw','0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin.nii.gz');
fgFileName   = fullfile(baseDir,'fibers','whole_brain_MRTRIX','right_occipital_MRTRIX.mat');

%% Initialize the Connectome
fe         = feConnectomeInit(dwiFile,dtFile,fgFileName);

%% Estimate the weights and isntall them in the fe structure
fefit      = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigfiber'),'sgdnn');
fe         = feSet(fe,'fit',fefit);

%% Now cross-validate the quality fo fit and install the result in the fe structure
fexval = feXvalidate(feGet(fe,'Mfiber'),feGet(fe,'dsigfiber'),'sgdnn');
fe         = feSet(fe,'xvalfit',fexval);

%% Now reduce the size of the fiber groups
fe = feConnectomeCull(fe);

%% Save it
feConnectomeSave(fe);

keyboard
return
