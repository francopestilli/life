function s_ms_fs_test_broca_wernicke_connection_FP(hemisphere)
%
% This script shows the basic workflow for a test of the conenctvity
% between two cortical areas.
%
% These are the steps involved in the process:
% (1) Load an ROI containing the two cortical areas (in this case a
%     putative Broca and putative Wernicke areas, were loaded from FS and 
%     combined into one ROI)
% (2) Load the vROI 
% (3) Load a whole-brain connectome.
% (4) Clip the whole-brain connectome within this vROI
% (5) Build the LIFE model for the clipped fiber group
% (6) Cull the connectome within the vROI
% (7) Perform an "and" operation between the fibers in the conenctome and
%     the cortical ROI.
% (8) If there are fibers left touching the cortical ROIs:
% (9) Find them and test the connection. This means test the increase in
%      RMSE in the vROI with and without the fascicles touching the 
%      cortical ROI.
%
% If there were no fibers left touchign the cortical ROI after stage (9).
% The test cannot be performed. An alternative could be that of ading more
% fibers touchign the cortical ROI to the connectome before (3). This can
% be performed by tracking-between the two corticla ROIs
%
%
% Written by Franco Pestilli (c) Stanford University 2013

% Paramters for the hemisphere, rois and white-matter volume
if notDefined('hemisphere'), hemisphere = 'rh';end
roiDir     = msPaths('bwrois');
bwRoiName  = fullfile(roiDir,sprintf('%s_bankssts_parsopercularis_parstriangularis',hemisphere));
bwRoiOperation = {'and both endpoints'};
vRoiName       = sprintf('%s_wm_roi_box_small',hemisphere);
saveDir = '~/Dropbox/';

% Parameters for the connectome and the tracking (these are used to decide
% which whole-brain connectome to load)
trackingType = 'prob';
lmax = 8;
bval = 2000;
rep  = 1;

% Paremters for saving/loading the FE structure
feSaveName = sprintf('%s_FE',vRoiName);
feSaveDir  = roiDir;

if ~exist(fullfile(feSaveDir,[feSaveName,'_culled.mat']),'file')
    if ~exist(fullfile(feSaveDir,[feSaveName,feSaveName,'.mat']),'file')
        % Load the Volume ROI
        vRoiFile   = fullfile(roiDir,vRoiName);
        vRoi       = dtiReadRoi(vRoiFile);
        
        % (3) Load the whole brain conectome
        [fgFileToLoad, ~] = msBuildFgFileName(trackingType,lmax,bval,rep);
        fg = fgRead(fgFileToLoad);
        
        % (4) Clip the whole-brain connectome to reside within the vROI
        maxVolDist   = 1; % The max distance in mm that we allow t node to be for accepting it into our vROI
        fgSaveName   = fullfile(roiDir,sprintf('%s_connectome',vRoiName));
        if ~exist(fgSaveName,'file')
            tic, fprintf('Clipping Connectome fibers that leave the volume ROI...\n');
            fg = feClipFibersToVolume(fg,vRoi.coords,maxVolDist);
            fgWrite(fg, fgSaveName,'pdb');
            fprintf('process completed in %2.3fhours\n',toc/60/60);
        end
        clear vRoi
        
        % (5) Build the life model for the clipped connectome
        dataRootPath  = fullfile('/biac2/wandell2/data/diffusion/pestilli/20120718_2975');
        subfolders    = fullfile('96dirs_b2000_1point5iso_1/');
        baseDir       = fullfile(dataRootPath,subfolders);
        dtFile        = fullfile(baseDir,'dt6.mat');
        dwiFile       = fullfile(dataRootPath,'preprocessed','run01_fliprot_aligned_trilin.nii.gz');
        dwiFileRepeat = fullfile(dataRootPath,'preprocessed','run02_fliprot_aligned_trilin.nii.gz');
        t1File        = fullfile(dataRootPath,'t1','t1.nii.gz');
        diffusionModelParams = [1,0];
        fe = feConnectomeInit(dwiFile,dtFile,fg,feSaveName,feSaveDir,dwiFileRepeat,t1File,diffusionModelParams);
        feConnectomeSave(fe,feSaveName)
        clear fg bwRoi
    else
        fprintf('[%s] Loading UNCULLED FE structure:\n%s\n',mfilename,fullfile(feSaveDir,[feSaveName,feSaveName,'.mat']))
        load(fullfile(feSaveDir,[feSaveName,feSaveName,'.mat']))    
    end
    
    % (6) Cull the connectome
    [fe, o] = feConnectomeCullNew(fe);
    save(fullfile(feSaveDir,[feSaveName,'_culled.mat']),'fe','o')
else % Load it 
    fprintf('[%s] Loading FE structure:\n%s\n',mfilename,fullfile(feSaveDir,[feSaveName,'_culled.mat']))
    load(fullfile(feSaveDir,[feSaveName,'_culled.mat']))
end

% (7) Perform an "and" operation between the fibers in the conenctome and
%     the cortical ROI.
fg = feGet(fe,'fibers acpc');
bwRoi      = dtiReadRoi(bwRoiName);
[fas, keepFascicles] = feSegmentFascicleFromConnectome(fg, {bwRoi}, bwRoiOperation, 'prob connectome');

% (8) Find them and test the connection. This means test the increase in
%     RMSE in the vROI with and without the fascicles touching the 
%     cortical ROI.
%
% Remove all the voxels from the connectome except the ones where the
% fascicle passes through. Fit the new model.
[feWithoutFas, feWithFas, connectivity] = feTestFascicle(fe,keepFascicles,0);

% Make a plot of the R-squared
WITH.r2       = median(feGetRep(feWithFas,   'vox  r2'));
WITH.rmse     = median(feGetRep(feWithFas,   'vox  rmse'));
WITH.rrmse    = median(feGetRep(feWithFas,   'vox  rmse ratio'));
WITH.rmseall  = (feGetRep(feWithFas,   'vox  rmse'));

WITHOUT.r2    = median(feGetRep(feWithoutFas,'vox  r2'));
WITHOUT.rmse  = median(feGetRep(feWithoutFas,'vox  rmse'));
WITHOUT.rrmse = median(feGetRep(feWithoutFas,'vox  rmse ratio'));
WITHOUT.rmseall  = (feGetRep(feWithoutFas,'vox  rmse'));

% Compute a test of the diference in rmse
% (a) Get the differece in rmse observed empiriclly
EmpiricalDiff = WITHOUT.rmse - WITH.rmse;

% (b) Compute the Null distribution by:
% (b.1) Combine all the rmse from both WITH and WITHOUT.
% (b.2) Compute 10,000 distributions of rmse one for WITH one for WITHOUT
%       with the same numerosity of the initial samples
% (b.3) Compute the difference between the medians of these 10,000
%       distributions.
NullSet     = [WITHOUT.rmseall WITH.rmseall];
sizeWith    = length(WITH.rmseall);
sizeWithout = length(WITHOUT.rmseall);
    
fprintf('[%s] Bootstrapping ...\n',mfilename)
nboots = 100000;
nullDistribution = nan(nboots,1);
parfor ibt = 1:nboots
    % Y = RANDSAMPLE(POPULATION,K)
    BootWith    = median(randsample(NullSet,sizeWith));
    BootWithout = median(randsample(NullSet,sizeWithout));
    nullDistribution(ibt) = BootWithout - BootWith;
end

% Plot the null distribution and the empirical difference
figName = sprintf('Test_Broca_Wernicke_rmse_%s_lmax%i',vRoiName,lmax);
fh = mrvNewGraphWin(figName);
[y,x] = hist(nullDistribution,100);
y = y./sum(y);
bar(x,y,'k')
hold on
plot([EmpiricalDiff,EmpiricalDiff],[0 max(y)],'r-','linewidth',2)
set(gca,'tickdir','out','box','off', 'ylim',[0 max(y)],'FontSize',16)
ylabel('Probability')
xlabel('Difference in rmse')

% (c) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
if max(nullDistribution)<EmpiricalDiff
    p = 100*1/nboots;
else
    p = sum(nullDistribution(sort(nullDistribution)>EmpiricalDiff));
end
dprime = EmpiricalDiff/std(nullDistribution);
title(sprintf('The probability of obtaining the difference by chance is less then %2.6f%%\nStrength of connectivity %2.3f',p,dprime), ...
    'FontSize',16)
saveFig(fh,fullfile(saveDir,figName))

% Make a scatter plto:
figName = sprintf('Test_Broca_Wernicke_rmse_SCATTER_%s_lmax%i',vRoiName,lmax);
fh = mrvNewGraphWin(figName);
hold on
plot([0 80],[0 80],'k-',[median(WITHOUT.rmse),median(WITHOUT.rmse)],[0 80],'k--', ...
     [0 80],[median(WITH.rmse),median(WITH.rmse)],'k--')
plot(WITHOUT.rmseall,WITH.rmseall,'o','color',[.3 .88 .7],'markerfacecolor',[.3 .88 .7]);
axis square
set(gca,'tickdir','out','box','off','xlim',[0 80],'ylim',[0 80],'FontSize',16)
ylabel('rmse WITH Connection')
xlabel('rmse WITHOUT Connection')
saveFig(fh,fullfile(saveDir,figName))

end % Main Function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)

printCommand = sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName);
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

% do the printing here:
eval(printCommand);
end

