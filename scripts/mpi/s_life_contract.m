function s_contrack
% function s_contrack
%
% This function uses the ctr*Batch*.m files to track fibers between V1 and
% LGN.
%
% It was a first attempt to track for the weiss fellowship.
%
% franco pestilli

dataRootPath = '/biac2/wandell6/data/arokem/ModelFits/';
subjDir      = 'FP20120420';
subfolders   = '150dirs_b1000_1';
dataDir      = fullfile(dataRootPath,subjDir,subfolders);

% initialize the parameters for the batch file.
ctrParams          = ctrInitBatchParams;

%% I. Set Naming and Directory Structure       
ctrParams.baseDir     = '/biac2/wandell6/data/arokem/WMRET/';                      
ctrParams.dtDir       = 'dti150trilin';
ctrParams.roiDir      = 'ROIs';

ctrParams.projectName = 'life_lgn_v1';
ctrParams.logName     = 'life_lgn_v1';
ctrParams.baseDir     = dataRootPath;
ctrParams.dtDir       = subfolders;
ctrParams.roiDir      = 'ROIs';

%% II. Set Subjects and ROIs
ctrParams.subs        = {subjDir};
ctrParams.roi1        = {'RLGN'};
ctrParams.roi2        = {'RV1_cleaned'};


%% III. Set Algorithm Parameters
ctrParams.nSamples     = 100000;    
ctrParams.maxNodes     = 240;      
ctrParams.minNodes     = 10;       
ctrParams.stepSize     = 1;        
ctrParams.pddpdfFlag   = 0;        
ctrParams.wmFlag       = 0;        
ctrParams.roi1SeedFlag = 'true';   
ctrParams.roi2SeedFlag = 'true';   
ctrParams.multiThread  = 0;
ctrParams.executeSh    = 0;

% use the parameters to created a batch contract file command.
[cmd infoFile] = ctrInitBatchTrack(ctrParams);

% track fibers...
system(cmd);

% create a batch command to score the fibers just created
numFibersToKeep = 1000; % returns the best fibers 
batchFileName  = ctrInitBatchScore(infoFile, numFibersToKeep, 1);

keyboard
% score the fibers
system(batchFileName);


