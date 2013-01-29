function [fe,feNew] = s_fe_delete_fibers
%
% Illustrate how to open up data and run a small linear fascicle evaluation
% (LIFE).
% 
% s_fe_delelte_fibers
%
% We should start making fePlot(), a gateway routine that takes an fe
% structure and some other arguments and generates visualizations of
% interesting quantities.
%
%  [uData, g] = fePlot(fe,plotType,varargin);
%
%
% Franco (c) 2012 Stanford VISTA team.


%% Initialize a small connectome for testing code
baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
dtFile     = fullfile(baseDir,'dti40','dt6.mat');
dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
fgFileName = fullfile(baseDir,'fibers','leftArcuateSmall.pdb');
fe         = feConnectomeInit(dwiFile,dtFile,fgFileName);

feFileName = fullfile(baseDir,'life',feGet(fe,'name'));
fefit      = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigfiber'),'sgdnn');
fe         = feSet(fe,'fit',fefit);
fexval     = feXvalidate(feGet(fe,'Mfiber'),feGet(fe,'dsigfiber'),'sgdnn');
fe         = feSet(fe,'xvalfit',fexval);

% This saves inside the dt6 dir
feConnectomeSave(fe,feFileName);

%% Delete some fibers:
% The select process builds the parameters inside of fe.life for the
% selected set of fibers.  The parameters inside of fe.life that are
% affected are
%   The M matrix
%   The fe.life.fibers entries
%
% We also store and keep the whole connectome in .fg.  But not all
% calculations with the M matrix use all the fibers.  In some cases we run
% smaller subsets.
%
nFibers      = feGet(fe,'n fibers');
fibersToKeep = 1:2:nFibers;
feNew        = feConnectomeSelectFibers(fe,fibersToKeep);

%% Estimate the weights
fefit = feFitModel(feGet(feNew,'Mfiber'),feGet(feNew,'dsigfiber'),'sgdnn');
feNew = feSet(fe,'fit',fefit);


%% Now cross-validate the quality fo fit
fexval = feXvalidate(feGet(feNew,'Mfiber'),feGet(feNew,'dsigfiber'),'sgdnn');
feNew = feSet(fe,'xvalfit',fexval);

keyboard
return
