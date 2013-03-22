%
% Illustrate how to open up data and run a small linear fascicle evaluation
% (LIFE).
%
% Written by Franco Pestilli (c) 2013 Stanford VISTA team.


% Select the type of fit done with or without L1 constrain
lType = 'L2'; % Option: L1 or L2

% Select the type of tractography
tractographyType = 'Prob'; % Tensor, Det or Prob

% Build a file name for the fe structure
feType = sprintf('fe%s%s',tractographyType,lType);

% Directory contaning the fe structures
fePath = '~/';

% Load the fe structure if not in memory
disp('Loading an FE structure')
if ~exist(feType,'var')
    load(fullfile(fePath,feType))
end

%% Read the fascicle
fg        = (feGet(eval(feType),'fibers acpc'));
allCoords = fefgGet(fg, 'uniqueimagecoords');

% Pick a voxel out of the ROI
vox = [6000:7000];% 800 1000
roiCoords = allCoords(vox,:);

%% Clip the connectome to be constrained within the volumeRoi.
tic, fprintf('Clipping Connectome fibers that leave the volume ROI...\n');
fg1 = feClipFibersToVolume(fg,roiCoords,.75);
fprintf('process completed in %2.3fhours\n',toc/60/60);

%% Initialize the Connectome
dwiFile    = feGet(eval(feType),'dwi file');
dwiFileRep = feGet(eval(feType),'dwi file rep');
dtFile     = feGet(eval(feType),'dt file');
t1file     = feGet(eval(feType),'anatomy file');
fe         = feConnectomeInit(dwiFile,dtFile,fg1,'test_singel_vox_fiber_removal',[],dwiFileRep,t1file);

%% Estimate the weights and install them in the fe structure
fe      = feSet(fe,'fit',feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn'));

% Now Remove half of the fascicles, the ones with highest weights
[w, keepFascicles] = sort(feGet(fe,'fiber weights'));
keepFascicles      = keepFascicles(floor(length(w)/1.6+1:end));
%feRedux           = feConnectomeReduceFibers(fe, zeros(size(keepFascicles)) );
feRedux            = feConnectomeReduceFibers(fe, ((keepFascicles)) );

% Now fit again
feRedux   = feSet(feRedux,'fit',feFitModel(feGet(feRedux,'Mfiber'),feGet(feRedux,'dsigdemeaned'),'sgdnn'));

% Compute some stats
r2Full     = 100*median(feGetRep(fe,'vox r2 '));
rmseFull   = median(feGetRep(fe,'vox rmse'));
rrmseFull  = median(feGetRep(fe,'vox rmse ratio'));
nW         = length(find(feGet(fe,'fiber weights'))) ;
Msiz       = size(feGet(fe,'mfiber'),1)/150;
ful = sprintf('\nFULL  connectome:  numFibers %i, numVox %i, ExpVar %2.3f%%, rmse %2.3f, Rrmse %2.3f\n', ...
        nW,Msiz,  r2Full,rmseFull,rrmseFull);
disp(ful)
    
r2Redux    = 100*median(feGetRep(feRedux,'vox r2 '));
rmseRedux  = median(feGetRep(feRedux,'vox rmse'));
rrmseRedux = median(feGetRep(feRedux,'vox rmse ratio'));
nWRedux    = length(find(feGet(feRedux,'fiber weights'))) ;
MsizRedux  = size(feGet(feRedux,'mfiber'),1)/150;
redux = sprintf('REDUX connectome:  numFibers %i, numVox %i,  ExpVar %2.3f%%, rmse %2.3f, Rrmse %2.3f\n', ...
        nWRedux,MsizRedux, r2Redux, rmseRedux,rrmseRedux);
disp(redux)

% Find the best fitting voxel:
[bestRmse,index] = sort(feGetRep(fe,'vox r2'));
rmseRedux        = (feGetRep(feRedux,'vox r2'));

% Compute the loss in r2
lossRmse = bestRmse - rmseRedux(index);

% Now sort the voxels by change in prediction given the removal of the
% fibers
[loss, lossIndex] = sort(lossRmse);

% Choose a voxel that had a large loss
% Plot the measured and predited signals
mrvNewGraphWin(sprintf('FULL %s model numVoxels %i',lType,length(vox)))
plot(feGetRep(fe,'dsig demeaned',index(4)'),'k.-')
hold on
plot(feGet(fe,'psigfiber',index(4)'),'r-')
set(gca,'box','off')
title(ful)

% Plot measured and predicted signal
mrvNewGraphWin(sprintf('REDUX %s model numVoxels %i',lType,length(vox)))
plot(feGetRep(feRedux,'dsig demeaned',index(4)'),'k.-')
hold on
plot(feGet(feRedux,'psigfiber',index(4)'),'r-')
title(redux)

% mrvNewGraphWin(sprintf('FULL vs. REDUX %s numVoxels %i',lType,length(vox)))
% plot(feGet(fe,'psigfiber'),feGet(feRedux,'psigfiber'),'ro')
% title('Correlation of predicted signals')


return



