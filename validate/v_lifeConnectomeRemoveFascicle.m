function fe = v_lifeConnectomeRemoveFascicle
%
% Generate a synthetic volume with the signal predicted by a Connectome.
% Illustrate how to use the fit of a Connectome to generate a prediction of signal in the DWI volume.
% 
% v_lifeConnectomeRemoveFascicle
%
% Franco (C) 2012 Stanford VISTA team.


%% Initialize a connectome
baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
dtFile     = fullfile(baseDir,'dti40','dt6.mat');
dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');
fgFileName = fullfile(baseDir,'fibers','leftArcuateSmall.pdb');

%% Initialize the Connectome
fe     = feConnectomeInit(dwiFile,dtFile,fgFileName);

%% Estimate the weights and install them in the fe structure
fefit  = feFitModel(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
fe     = feSet(fe,'fit',fefit);

%% Now cross-validate the quality fo fit and install the result in the fe structure
fexval = feXvalidate(feGet(fe,'Mfiber'),feGet(fe,'dsigdemeaned'),'sgdnn');
fe     = feSet(fe,'xvalfit',fexval);

%% Now selected a subset of fibers to romove from the connectome.
% Identify the indices of the fibers to remove.
fibersToRemove.percent = 10;
fibersToRemove.num     = round(feGet(fe,'nfibers') .* (fibersToRemove.percent/100));
fibersToRemove.indices = sort(randsample(feGet(fe,'nfibers'),fibersToRemove.num));

%% Predict the signal WITH the fibers
rmse_with    = (feGet(fe, 'vox rmse'));

%% Predict the signal WITHOUT the fibers
rmse_without    = (feGet(fe,     'vox rmse test',fibersToRemove.indices));

%% make a scatter plot of the loss in explained variance.
mrvNewGraphWin('Increase in RMSE without optic radiation fibers')
plot(rmse_with,rmse_without,'ro','MarkerFaceColor','r')
axlim = get(gca,'xLim');axis square;set(gca,'yLim',axlim);
hold on, 
plot([axlim(1);axlim(2)], [axlim(1); axlim(2)],'k-')
ylabel('RMSE without some fibers');xlabel('RMSE with all fibers')

%% make a histogram RMSE
mrvNewGraphWin('RMSE with and without optic radiation fibers')
[y,x] = hist(rmse_without,200);bar(x,y,'r','BarWidth',0.65)
hold on
[y,x] = hist(rmse_with,200);bar(x,y,'k','BarWidth',0.35)
ylabel('Occurrences');xlabel('RMSE');
legend({'Without some fibers','With all fibers'})

%% Generate maps of RMSE 
niftiName = fullfile(feGet(fe,'savedir'),'rmse_with_new_remove_test');
rmse_w_map = feValues2volume(rmse_with,feGet(fe,'roi coords'),feGet(fe,'map size'));
nii         = niftiCreate('data',rmse_w_map,...
                  'qto_xyz',feGet(fe,'xform img 2 acpc'),...
                  'fname',niftiName);
niftiWrite(nii);

niftiName = fullfile(feGet(fe,'savedir'),'rmse_without_new_remove_test');
rmse_w_map = feValues2volume(rmse_without,feGet(fe,'roi coords'),feGet(fe,'map size'));
nii         = niftiCreate('data',rmse_w_map,...
                  'qto_xyz',feGet(fe,'xform img 2 acpc'),...
                  'fname',niftiName);
niftiWrite(nii);

return


%%
% Render a nice image of the fibers using build3DFrame and
% buildSurfaceFromFrame
fePlot(fe,'connectome');



