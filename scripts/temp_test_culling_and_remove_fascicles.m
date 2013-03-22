if notDefined('trackingType'),trackingType = 'tensor';end
if notDefined('lmax'),        lmax         = 6;end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = 1;end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_figures');end

% High-resolution Anatomy
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';

% Path to the ROIs used to segment the fascicle
roisPath = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/ROIs_restric_connectomes';

% These are the ROIs used to segment the fascicle
rois = {fullfile(roisPath,'ROI_Anterior.mat'),  ...
        fullfile(roisPath,'ROI_Inferior.mat'),  ...
        fullfile(roisPath,'ROI_posterior.mat'),  ...
        fullfile(roisPath,'ROI_Superior.mat')};

% These are the logical operations to apply to each ROI
operation = {'and','not','and','not'};

% Information on the path to the files to load.
% This is where the inputs will be loaded from
feFileToLoad = msBuildFeFileName(trackingType,lmax,bval,rep,diffusionModelParams);

% Get the fe structure
disp('loading the LiFE structure...')
% load(feFileToLoad);
% 
% [feTensorL1, oTensorL1] = feConnectomeCull(fe,1000,'sgdl1nn');
% [feTensorL2, oTensorL2] = feConnectomeCull(fe,1000,'sgdnn');

%save ~/feTensor.mat feTensorL1 feTensorL2 oTensorL1 oTensorL2 

trackingType = 'probabilistic';

% Information on the path to the files to load.
% This is where the inputs will be loaded from
feFileToLoad = msBuildFeFileName(trackingType,lmax,bval,rep,diffusionModelParams);

% Get the fe structure
disp('loading the LiFE structure...')
load(feFileToLoad);

[feProbL2, oProbL2] = feConnectomeCull(fe,1000,'sgdnn');
[feProbL1, oProbL1] = feConnectomeCull(fe,1000,'sgdl1nn');

save ~/feProb.mat feProbL1 feProbL2 oProbL1 oProbL2 

trackingType = 'deterministic';

% Information on the path to the files to load.
% This is where the inputs will be loaded from
feFileToLoad = msBuildFeFileName(trackingType,lmax,bval,rep,diffusionModelParams);

% Get the fe structure
disp('loading the LiFE structure...')
load(feFileToLoad);

[feDetL2, oDetL2] = feConnectomeCull(fe,1000,'sgdnn');
[feDetL1, oDetL1] = feConnectomeCull(fe,1000,'sgdl1nn');

save ~/feDet.mat feDetL1 feDetL2 oDetL1 oDetL2 
keyboard


% % Extract the connectome
% fgP = feGet(feP,'fibers acpc');
% feConnectomeDisplay( feSplitLoopFibers(fgP), figure ) 
% 
% % Perform a series of end and not operations to segement the OR/ILF
% [fgPF keepFG] = feSegmentFascicleFromConnectome(fgP, rois, operation, 'fgNF');
% 
% % Now remode the fibers of the fascicle from the LiFE model.
% feNotPF = feConnectomeReduceFibers(feP,~keepFG);
% 
% % Now re-fit the model.
% feNotPF = feSet(feNotPF,'fit',feFitModel(feGet(feNotPF,'M fiber'),feGet(feNotPF,'dsig demeaned'),'sgdnn'));
% feConnectomeDisplay( feSplitLoopFibers(feGet(feNotPF,'fibers img')), figure ) 
% 
% % Now remode the fibers of the fascicle from the LiFE model.
% fePF = feConnectomeReduceFibers(feP,keepFG);
% 
% % Now re-fit the model.
% fePF = feSet(fePF,'fit',feFitModel(feGet(fePF,'M fiber'),feGet(fePF,'dsig demeaned'),'sgdnn'));
% feConnectomeDisplay( feSplitLoopFibers(feGet(fePF,'fibers img')), figure ) 
% 
% % Test if combining the two fiber groups we win anything.
% fgBoth = fgMerge(fgP,fgN,'combined tensor and lmax12 prob');

% Re build the LiFE structure with only the good fibers
% dwiFile         = feGet(fe,'dwifile') ;
% dtFile          = feGet(fe,'dtfile');
% feFileName      = [feGet(fe,'name'),'bestFG'];
% savedir         = feGet(fe,'savedir');
% dwiFileRepeated = feGetRep(fe,'dwifile');
% anatomyFile     = feGet(fe,'anatomy file'); 
% tensorModel     = feGet(fe,'model tensor');
% 
% clear fe*
% clear fgPF fgN
% 
% fe = feConnectomeInit(dwiFile,dtFile,fgC,feFileName,savedir,dwiFileRepeated,anatomyFile,tensorModel(1:2));
% 
% feBoth = feConnectomeCull(feBoth);
% 
% r = median(feGetRep(feBoth,'vox rmse ratio'));
% r2 = feGetRep(fen,'total r2');
% save ~/Dropbox/rBoth.mat r r2
% 
keyboard