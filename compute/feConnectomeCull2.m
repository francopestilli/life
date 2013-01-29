function [fe fefit fexval] = feConnectomeCull2(fe,minWeight, prctR2redux,maxNumIter)
%
% This function si currently not being used. But it might be useful in the
% future releases.
%
% Find the smallest set fibers in a connectome that explain most of the variance.
%
% It iteratively finds and drops fibers with smallest weights (zero weights).
%
%   [fe fefit fexval] = feConnectomeCull2(fe,maxNumInter, minWeight, prctR2redux,maxNumIter)
% 
% Inputs:
%   fe            - An fe structure, see feCreate.m or v_lifeExample.m
%   minWeight     - The lowest weight accepted.  Default, is zero.
%   prctR2redux   - The max reduction in percent variance explained that we
%                   allow before exiting. Default is zero.
%   maxNumIter    - Maximum number of iterations allowed before retuning
%                   the current connectome. Default 1000;
%
% Outputs: 
%   fe            - The fe structure with  the culled fibers.
%
% Example:
%   baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
%   dtFile     = fullfile(baseDir,'dti40','dt6.mat');
%   dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
%   fgFileName = fullfile(baseDir,'fibers','leftArcuateSmall.pdb');
%   fe         = feConnectomeInit(dwiFile,dtFile,fgFileName);
%   fe         = feConnectomeCull(fe);
% 
% Franco (c) 2012 Stanford VISTA Team.

% This is the percent of fibers that are removed at each iteration. The
% larger this number the faster the convergence, but the more likely to
% delete important fibers and produce a reduced connectome with a larger
% loss of R2.
if notDefined('minWeight'), minWeight     = 0;end

% The percent change in R2 from that of the orignal model. 
% The smaller this number the faster the convergence, the more the fibers
% kept in the connectome.
if notDefined('prctR2redux'), prctR2redux = 0;end

% Maximum number of iterations.
if notDefined('maxNumIter'), maxNumIter   = 1000;end

% Show a plot of the results:
makePlot = 0;

% Bookkeeping variables.
all_orig_r2 = [];
all_curr_r2 = [];
all_nfibers = [];

% Fit the model and then reduced it by only accepting fibers that pass
% the minWeights threshold.
for iter = 1:maxNumIter
  fefit        = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dSig demeaned'),'sgdnn');  
   
  % Check whether we start loosing percent variance explained when
  % crossvalidating
  fexval = feXvalidate(feGet(fe,'Mfiber'),feGet(fe,'dSig demeaned'),'sgdnn');
  
  % STOP criteria: When the fit quality decreases we stop culling.
  % Keep track of the quality of fit.
  if iter == 1
    originalR2 = fexval.r2.mean;  
    currentR2  = fexval.r2.mean;
  
    % Book-keeping
    all_orig_r2 = [all_orig_r2, originalR2];
    all_curr_r2 = [all_curr_r2,  currentR2];
    all_nfibers = [all_nfibers, feGet(fe,'n fibers')];
  else
    currentR2  = fexval.r2.mean;
    % Check whether the new R2 is higher then the original. If it is, set
    % the target R2 to the current one.
    if currentR2 > originalR2, originalR2 = currentR2;end
  end
  
  if currentR2 < ((originalR2) - (originalR2*(prctR2redux/100)))
    % Book-keeping
    all_orig_r2 = [all_orig_r2, originalR2];
    all_curr_r2 = [all_curr_r2,  currentR2];
    all_nfibers = [all_nfibers, length(find(fefit.weights > minWeight))];    
    break

  end
  fprintf('\n\n[%s] n iter: %i, Original R^2: %2.3f, Current R^2: %2.3f.\n\n',mfilename, iter,originalR2,currentR2)
  
  % Compute the cut-off for the fibers
  % minWeight    = prctile(fefit.weights,minWeight);
  fibersToKeep = find(fefit.weights > minWeight);
  
  % Reduce the connectome.
  fe = feConnectomeSelectFibers(fe,fibersToKeep);
  
  % Book-keeping
  all_orig_r2 = [all_orig_r2, originalR2];
  all_curr_r2 = [all_curr_r2,  currentR2];
  all_nfibers = [all_nfibers, feGet(fe,'n fibers')];  
end

fprintf('[%s] Done culling, in %i iterations.\n    Original R^2: %2.3f, Current R^2: %2.3f.\n',mfilename, iter,originalR2,currentR2)

if makePlot
  % Make a plot of the culling results
  mrvNewGraphWin('Culling procedure');
  subplot(1,2,1),plot(all_curr_r2,'ko-','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',16)
  hold on
  subplot(1,2,1),plot(all_orig_r2,'ro-','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',12)
  subplot(1,2,1),plot([iter iter],get(gca,'ylim'))
  ylabel('Percent Variance Explaied (R^2)'), xlabel('Interation'),legend({'Current R^2', 'Target R2'})
  subplot(1,2,2),plot(all_nfibers,'k.-');
  set(gca,'yLim',[min(all_nfibers)-5 max(all_nfibers) + 5]), ylabel('Number of Fibers'), xlabel('Iteration');
  hold on
  subplot(1,2,2),plot([iter iter],get(gca,'ylim'))
end

return