function s_life_fg_afq_stt_vs_tend_vs_mrtrix(runType,fitType,save)
% function s_life_fg_afq_stt_vs_tend_vs_mrtrix
%
% Loads some tracks produced with STT or TEND. 
% 
% Then produces the quality of fit for the tracks, before and after fitting with LiFE.
%
% Can run the fits on the same ROI (voxels common to the two fibers
% geenrated by the two algorithms), or on the independent unique voxels
% identified by each algorithm's fiber groups.
%
% Franco (C) 2012 Stanford VISTA team.

if notDefined('dtiDataType'); dtiDataType = 'b1000'; end
if notDefined('lambda');           lambda = 0;       end
if notDefined('fitType');         fitType = 'sgd';   end
if notDefined('runType');         runType = 'i';     end % Run life on common voxels defined as the intersection between the two fibers. 
if notDefined('save');               save = 1;       end % Save or not files and figures
if notDefined('which_fascicles');which_fascicles = 1:20;end % There are 20 fascicles

baseDir  = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/';
dtDir    = fullfile('150dirs_b1000_1');
dtFile   = fullfile(baseDir,dtDir,'dt6.mat');
dwiFile  = fullfile(baseDir,'raw',sprintf('0009_01_DWI_2mm150dir_2x_%s_aligned_trilin.nii.gz',dtiDataType));

fiberDir = fullfile(baseDir,dtDir,'fibers','whole_brain');
saveDir  = fullfile(baseDir,'LiFE','compare_stt_tend_mrtrix');
saveName = sprintf('%s_%s_%s_%s',mfilename,dtiDataType,fitType,runType);

% Set default parameters
[stt tend mrtrix nifti dwi] = mctSetDefaultParms(dtiDataType,fiberDir,dtFile,dwiFile,baseDir,fitType,lambda,saveDir,saveName,save);

for fg_i = 1:length(which_fascicles)
  fprintf('\n[%s] Running tests on fascicle #%i...\n',mfilename,fg_i);
  
  % Check how we want to run LiFE
  switch runType
    case {'independently', 'i','ind'} % On independent voxels for TEND and STT
      disp('Running LiFE on independent voxels between two algorithms...')
      stt.runType  = runType;
      tend.runType = runType; 
      mrtrix.runType =runType;
      [stt tend mrtrix] = runLifeBothAlgos(tend,stt,mrtrix,dwi,nifti,fg_i);
      
    case {'same','common','s','c'}    % on the same, common voxels between the two fiber groups
      disp('Running LiFE on common voxels between two algorithms...')
      stt.runType  = runType;
      tend.runType = runType;
      mrtrix.runType =runType;
      [stt tend mrtrix] = runLifeCommonVoxels(tend,stt,mrtrix,dwi,nifti,fg_i);
      
    otherwise
      keyboard
  end
  
  drawnow; close all; drawnow
  fprintf('\n[%s] DONE running tests on fascicle #%i...\n',mfilename,fg_i);
end

end


%%%%%%%%%%%%%%%%%%%%%%
% mctSetDefaultParms %
%%%%%%%%%%%%%%%%%%%%%%
function [stt tend mrtrix nifti dwi] = mctSetDefaultParms(dtiDataType,fiberDir,dtFile,dwiFile,baseDir,fitType,lambda,saveDir,saveName, save)
%
% Set the dafault parameters for the two algorithms being tested.
%
% Franco

% load the dwi and dti files        
dwi          = dwiLoad(dwiFile);
[dtiF, dtiH] = mrDiffusion('off',dtFile);

% Set parameters for the STT fiber group
stt.algo           = 'STT';
stt.dtiDataType    = dtiDataType;
stt.fg.name        = fullfile(sprintf('%s_%s',fiberDir,stt.algo),sprintf('MoriGroups_clean_D5_L4_%s.mat',stt.algo));
stt.xform.img2acpc = dtiGet(dtiH,'img 2 acpc xform');
stt.xform.acpc2img = dtiGet(dtiH,'acpc 2 img xform');
stt.fitType        = fitType;
stt.lambda         = lambda;
stt.savedir        = [saveDir,'_',stt.algo];
stt.save           = save; % save figures and results
stt.savename       = saveName;

% Set parameters for the TEND fiber group
tend.algo           = 'TEND';
tend.dtiDataType    = dtiDataType;
tend.fg.name        = fullfile(sprintf('%s_%s',fiberDir,tend.algo),sprintf('MoriGroups_clean_D5_L4_%s.mat',tend.algo));
tend.xform.img2acpc = dtiGet(dtiH,'img 2 acpc xform');
tend.xform.acpc2img = dtiGet(dtiH,'acpc 2 img xform');
tend.fitType        = fitType;
tend.lambda         = lambda;
tend.savedir        = [saveDir,'_',tend.algo];
tend.save           = save; % save figures and results
tend.savename       = saveName;

% Set parameters for the TEND fiber group
mrtrix.algo           = 'MRTRIX';
mrtrix.dtiDataType    = dtiDataType;
mrtrix.fg.name        = fullfile(sprintf('%s_%s',fiberDir,mrtrix.algo),sprintf('MoriGroups_clean_D4_L3_%s.mat',mrtrix.algo));
mrtrix.xform.img2acpc = dtiGet(dtiH,'img 2 acpc xform');
mrtrix.xform.acpc2img = dtiGet(dtiH,'acpc 2 img xform');
mrtrix.fitType        = fitType;
mrtrix.lambda         = lambda;
mrtrix.savedir        = [saveDir,'_',mrtrix.algo];
mrtrix.save           = save; % save figures and results
mrtrix.savename       = saveName;

% Close mrDiffusion
close(dtiF), drawnow

% load the T1, this is used to display a slice with the fiber groups
disp('Loading T1')
nifti = readFileNifti(fullfile(baseDir,'t1','t1.nii.gz'));

end


%%%%%%%%%%%%%%%%%%%%%%%
% runLifeCommonVoxels %
%%%%%%%%%%%%%%%%%%%%%%%
function [stt tend mrtrix] = runLifeCommonVoxels(tend,stt,mrtrix,dwi,nifti,fg_i)
%
% Runs life on the common voxels between the stt and tend fibers.
%
% Franco

% Load the STT or TEND fibers and find the common voxels.
fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,tend.fg.name) 
temp_fg      = dtiLoadFiberGroup(tend.fg.name); 
tend.fg.acpc = temp_fg(fg_i);clear temp_fg;
tend.fg.img  = dtiXformFiberCoords(tend.fg.acpc, tend.xform.acpc2img,'img');
tend.coords  = fefgGet(tend.fg.img, 'uniqueimagecoords');

fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,stt.fg.name)
temp_fg      = dtiLoadFiberGroup(stt.fg.name); 
stt.fg.acpc = temp_fg(fg_i);clear temp_fg;
stt.fg.img  = dtiXformFiberCoords(stt.fg.acpc, stt.xform.acpc2img,'img');
stt.coords  = fefgGet(stt.fg.img, 'uniqueimagecoords');

fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,mrtrix.fg.name)
temp_fg      = dtiLoadFiberGroup(mrtrix.fg.name); 
mrtrix.fg.acpc = temp_fg(fg_i);clear temp_fg;
mrtrix.fg.img  = dtiXformFiberCoords(mrtrix.fg.acpc, mrtrix.xform.acpc2img,'img');
mrtrix.coords  = fefgGet(mrtrix.fg.img, 'uniqueimagecoords');

% find the unique common coordinates between the two fiber groups.
coords = intersect(stt.coords,tend.coords,'rows');
coords = intersect(coords,mrtrix.coords,'rows');

% Now set the coordinates for both algorithms to the common ones.
stt.coords  = coords;
tend.coords = coords;
mrtrix.coords = coords;

% Now run life for the two algorithsm.
[stt tend mrtrix] = runLifeBothAlgos(tend,stt,mrtrix,dwi,nifti,fg_i);

end


%%%%%%%%%%%%%%%%%%%%
% runLifeBothAlgos %
%%%%%%%%%%%%%%%%%%%%
function [stt tend mrtrix] = runLifeBothAlgos(tend,stt,mrtrix,dwi,nifti,fg_i)
%
% This function runs life for the two algorithms being tested.
%
% Franco

% Run LiFE
tend = mctRunLifeFitsLocal(tend,dwi,fg_i);

% Plot results and fiber groups
tend = mctPlotResultsLocal(tend,nifti,fg_i);

% Run LiFE
stt = mctRunLifeFitsLocal(stt,dwi,fg_i);

% Plot results and fiber groups
stt = mctPlotResultsLocal(stt,nifti,fg_i);

% Run LiFE
mrtrix = mctRunLifeFitsLocal(mrtrix,dwi,fg_i);

% Plot results and fiber groups
mrtrix = mctPlotResultsLocal(mrtrix,nifti,fg_i);
end

  
%%%%%%%%%%%%%%%%%%%%
% runLifeFitsLocal %
%%%%%%%%%%%%%%%%%%%%
function val = mctRunLifeFitsLocal(val,dwi,fg_i)
%
% Performs a run of the life model and evaluates it against the
% original model
%
% Franco

% Load the contrack fibers and the STT or TEND fibers and add them together.
fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,val.fg.name)
temp_fg      = dtiLoadFiberGroup(val.fg.name); 
val.fg.acpc = temp_fg(fg_i);clear temp_fg;
val.fg.img  = dtiXformFiberCoords(val.fg.acpc, val.xform.acpc2img,'img');

% Only compute the ROI coordinates for the independent ROI case.
switch val.runType
  case {'independently', 'i','ind'}
    val.coords  = fefgGet(val.fg.img, 'uniqueimagecoords');
  otherwise
end
fprintf('\n[%s] the ROI used has size(%i,%i)\n',mfilename,size(val.coords))

% Build the MicroTrack model
[Afiber, Aiso, val.dSig, dSig_demeaned] = mctBuildDiffusionModel(dwi,val.fg.img,val.coords,[],'ones');
val.nFiber = size(Afiber,2); 

% build the full model
Afull = [Afiber, Aiso];

% Fit the model.
val.life.w.fiber = mctFitDiffusionModel(Afiber, dSig_demeaned, val.fitType,val.lambda);

% Fit the isotropic
val.life.w.iso  = Aiso \ val.dSig;
val.life.w.full = [val.life.w.fiber; val.life.w.iso];

% Predict the signal 
val.life.pSig  = mctComputePredictedSignal(Afull,val.life.w.full);

% Overall quality across voxels
[val.life.rmse, val.life.r2] = mctComputePredictionQuality(val.dSig, val.life.pSig,1);

% Compute residuals, this is done on the demeaned (fiber) signal. 
% This signal was not explained by the model 
val.life.resSig = (val.dSig - val.life.pSig) + (val.dSig - dSig_demeaned);

% Compute the diffusion signal predicted by the original fiber model.
val.nVoxels = size(Aiso,2); 
val.nBvecs  = size(Afiber,1) / val.nVoxels;
Aorig       = [sum(Afiber,2)  Aiso * val.life.w.iso]; 
w           = Aorig \ val.dSig;
val.orig.w.fiber = w(1);
val.orig.w.iso   = w(2);
val.orig.w.full  = w;
val.orig.pSig    = mctComputePredictedSignal(Aorig,val.orig.w.full);
[val.orig.rmse, val.orig.r2] = mctComputePredictionQuality(val.dSig, val.orig.pSig,1);

val.orig.resSig = (val.dSig - val.orig.pSig) + Aorig(:,2);

% Select the fibers explaining most of the varicance
val.life.w.criterion(1) = 0.75; % get the top 75% fibers
val.life.w.criterion(2) = quantile(val.life.w.fiber,val.life.w.criterion(1));
val.life.hits        = val.life.w.fiber >= val.life.w.criterion(2); % the fibers with large weights
val.life.fa          = ~val.life.hits;

%fprintf('\n[%s] Extracting good and bad fibers...\n',mfilename)
%val.fg.good = fgExtract(val.fg.acpc,find(val.life.hits),'keep');
%val.fg.bad  = fgExtract(val.fg.acpc,find(val.life.fa),'keep');

end


%%%%%%%%%%%%%%%%%%%%%%%
% mctPlotResultsLocal %
%%%%%%%%%%%%%%%%%%%%%%%
function val = mctPlotResultsLocal(val,nifti,fg_i)
%
% plots results of a life fit.
%
% Franco

% plotting fiber predictions against the signal  
val.figs.hFit = mctDisplayModelFit(val.life.w.fiber,val.life.w.iso, val.dSig, val.life.pSig, val.life.r2,  ...
                                                                              val.orig.pSig, val.orig.r2, val.fitType);
% plot quality of fit.
val.figs.hQual = mctPlotQualityOfFit(val.life.r2,val.orig.r2, val.life.pSig,val.orig.pSig, val.dSig);

% Plot fibers with a slice.
%val.figs.hGood(1) = mctNfgDisplayStrands(val.fg.good,[.3 .2 .75],[],[],[],20,.1.*ones(size(val.fg.good.fibers)));
%hold on
%val.figs.hGood(2) = mctDisplayBrainSlice(nifti,[0 0 -7]);
%title('Good fibers')

%val.figs.hBad(1) = mctNfgDisplayStrands(val.fg.bad,[.75 .2 .3],[],[],[],20,.1.*ones(size(val.fg.bad.fibers)));
%hold on
%val.figs.hBad(2) = mctDisplayBrainSlice(nifti,[0 0 -7]);
%title('Bad fibers')

%val.figs.hAll(1) = mctNfgDisplayStrands(val.fg.acpc,[.2 .75 .3],[],[],[],20, 8 .* val.life.w.fiber);
%hold on
%val.figs.hAll(2) = mctDisplayBrainSlice(nifti,[0 0 -7]);
%title('All fibers scaled by weight')

if val.save
  val = feSaveResultsLocal(val);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mctSaveLifeResultsLocal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = feSaveResultsLocal(val)
%
% Save the results to disk.
%
% Franco

% Remove blank spaces in files name.
val.fg.acpc.name(isspace(val.fg.acpc.name))='_';
val.file.name  = sprintf('%s_%s_%s',val.savename,val.algo,val.fg.acpc.name);
val.file.results = fullfile(val.savedir,val.file.name);
fprintf('\n\nSaving file: %s\n\n',val.file.results)
save(val.file.results,'-v7.3')

disp('Saving results figures ... ')
% eps for the plots
% saveas(val.figs.hFit ,fullfile(val.savedir,[val.file.name,'_Fit_',val.algo]),'fig')
savefigvista(val.figs.hFit ,[val.file.name,'_Fit_',val.algo],'eps',val.savedir,1,0);

% saveas(val.figs.hQual ,fullfile(val.savedir,[val.file.name,'_Qual_',val.algo]),'fig')
savefigvista(val.figs.hQual ,[val.file.name,'_Qual_',val.algo],'eps',val.savedir,1,0);

if isfield(val.figs,'hAll')
  
  % fig for the fibers
  disp('Saving fibers figures ... ')
  saveas(val.figs.hAll(1) ,fullfile( val.savedir,[val.file.name,'_FgAll_',val.algo]),'fig')
  saveas(val.figs.hGood(1) ,fullfile(val.savedir,[val.file.name,'_FgGood_',val.algo]),'fig')
  saveas(val.figs.hBad(1) ,fullfile( val.savedir,[val.file.name,'_FgBad']),'fig')
end

end

