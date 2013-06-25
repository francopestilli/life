function s_ms_mt_ips()
%
% This script performs a test of conenctivity of MT+ (LO1 and LO2) with
% IPS0. THe following are the steps we perform:
%  - It loads a culled half hemisphere connectome (FE structure). 
%  - It loads the MT+ ROI.
%  - It loads the IPS0 ROI. 
%  - It finds the fascicles connecting MT+ and IPS0
%  - It reduces the connectome to the voxels and fibers of the conenctions
%  - It Perform a bootstrap test WITH and WITHOUT the connection between
%    MT+ and IPS0.
%
% Written by Franco Pestilli (c) Stanford University, Vista Team 2013
% High-resolution Anatomy
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';
t1     = niftiRead(t1File);

% Load the FE structure
feDir      = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/fe_hemispheres';
feFileName = 'fe_culled_FP_150_B2000_LMAX8_left.mat';
fprintf('[%s] Loading the FE structure...\n',filename)
load(fullfile(feDir,feFileName));

% Load the ROIs
roiDir = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/mt_roi';
mtFileName   = fullfile(roiDir,'LTO1_O2_cleaned.mat');
ips0FileName = fullfile(roiDir,'LIPS0_cleaned.mat');
mt   = dtiReadRoi(mtFileName);
ips0 = dtiReadRoi(ips0FileName);

% Find the fascicles in the connectome that touch both MT+ and IPS0.
operations = {'and endpoints','and endpoints'};
rois       = {mt,ips0};
fg         = feGet(fe,'fg acpc');
[mtIpsFG, keepFascicles] = feSegmentFascicleFromConnectome(fg, rois, operations, 'mt_ips_zero');
%[figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(mtIpsFG.fibers);

% Clean the fibers by length, fibers that too long are likely to go far
% frontal and not just touch MT+ and IPS0.
[Lnorm, Lmm] = mbaComputeFiberLengthDistribution(mtIpsFG, true);
maxSD        = 1; % Max standard deviation of the fibers to keep in the group.
fibers2delete  = Lnorm > maxSD;

% Now let's get the indices of the fibers in the FE structure:
fasIndices = find(keepFascicles);
fasIndices = fasIndices(fibers2delete);

% Now let's mark the fascicles as deleted.
keepFascicles(fasIndices) = false;

% Show te new fiber group
%mtIpsFG = fgExtract(mtIpsFG,Lnorm<maxSD,'keep');
%h = mbaDisplayBrainSlice(t1, [-18 0 0]);
%hold on
%h = mbaDisplayBrainSlice(t1, [0 0 -14 ]);
%h = mbaDisplayBrainSlice(t1, [0 -40 0]);
%[figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(mtIpsFG.fibers,gcf);
%delete(lightHandle)
%lightHandle = camlight('left');
%view(-75,30);
%axis([-67 -18 -110 -40 -18 80]);

% Reduce the connectome to the voxels and fibers of the connection between
% MT+ and IPS0.
[feWithoutFas, feWithFas, ~] = feTestFascicle(fe,keepFascicles,0);

% Make a plot of the R-squared
WITH.rmse      = median(feGetRep(feWithFas,   'vox  rmse'));
WITH.rrmse     = median(feGetRep(feWithFas,   'vox  rmse ratio'));
WITH.rmseall   = (feGetRep(feWithFas,   'vox  rmse'));

WITHOUT.rmse    = median(feGetRep(feWithoutFas,'vox  rmse'));
WITHOUT.rrmse   = median(feGetRep(feWithoutFas,'vox  rmse ratio'));
WITHOUT.rmseall = (feGetRep(feWithoutFas,'vox  rmse'));

% Compute a test of the diference in rmse
% (1) Get the differece in rmse observed empiriclly
EmpiricalDiff = WITHOUT.rmse(1) - WITH.rmse(1);

% (2) Compute the Null distribution by:
% (2.1) Combine all the rmse from both WITH and WITHOUT.
% (2.2) Compute 10,000 distributions of rmse one for WITH one for WITHOUT
%       with the same numerosity of the initial samples
% (2.3) Compute the difference between the medians of these 10,000
%       distributions.
NullSet     = [WITHOUT.rmseall WITH.rmseall];
sizeWith    = length(WITH.rmseall);
sizeWithout = length(WITHOUT.rmseall);

nboots = 20000;
nullDistribution = nan(nboots,1);
parfor ibt = 1:nboots
    % Y = RANDSAMPLE(POPULATION,K)
    BootWith    = median(randsample(NullSet,sizeWith));
    BootWithout = median(randsample(NullSet,sizeWithout));   
    nullDistribution(ibt) = BootWithout - BootWith;
end

% Plot the null distribution and the empirical difference
figName = sprintf('Test_MT_IPS0_connection_rmse_%s',feFileName(1:end-4));
fh = mrvNewGraphWin(figName);
[y,x] = hist(nullDistribution,100);
y = y./sum(y);
bar(x,y,'k')
hold on
plot([EmpiricalDiff,EmpiricalDiff],[0 max(y)],'r-','linewidth',2)
set(gca,'tickdir','out','box','off', 'ylim',[0 max(y)],'FontSize',16)
ylabel('Probability')
xlabel('Difference in rmse')

% (3) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
if max(nullDistribution)<EmpiricalDiff
    p = 100*1/nboots;
else
    p = sum(nullDistribution(sort(nullDistribution)>EmpiricalDiff));
end
title(sprintf('The probability of obtaining the difference by chance is less then %2.3f%%',p), ...
    'FontSize',16)
saveFig(fh,fullfile(saveDir,figName),1)
 
% Make a scatter plot:
figName = sprintf('Test_MT_IPS0_connection_rmse_SCATTER_%s',feFileName(1:end-4));
fh = mrvNewGraphWin(figName);
hold on
plot([0 80],[0 80],'k-',[median(WITH.rmse),median(WITH.rmse)],[0 80],'k--',[0 80],[median(WITHOUT.rmse),median(WITHOUT.rmse)],'k--')
plot(WITH.rmseall,WITHOUT.rmseall,'o','color',[.86 .65 .55],'markerfacecolor',[.86 .65 .55]);
axis square
set(gca,'tickdir','out','box','off','xlim',[0 80],'ylim',[0 80],'FontSize',16)
xlabel('rmse WITH Tract')
ylabel('rmse WITHOUT Tract')
saveFig(fh,fullfile(saveDir,figName),1)

% Make a plot of the effect size
figName = sprintf('Test_MT_IPS0_connection_rmse_SCATTER_%s',feFileName(1:end-4));
fh = mrvNewGraphWin(figName);
hold on

