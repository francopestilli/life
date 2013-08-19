function s_ms_mt_ips_tracts_figures(hemisphere,saveDir)
%
% This script shows the tracts conencting to MT+ before and after cullng each hemisphere.
% Also it shows the tracts conencting MT+ (LO1 and LO2) with
% IPS0. 
%
% Written by Franco Pestilli (c) Stanford University, Vista Team 2013

if notDefined('saveDir'), 
    saveDir = '~/Dropbox/connectomes_plot_mt_ips/';
end
if notDefined('hemisphere'), hemisphere = 'left';end
tractColor = [.95 .9 .86];%[.3 .7 .9];

% Set all the variables that depend on the hemisphere
switch hemisphere
    case {'left'}
        feFileName = 'fe_culled_FP_150_B2000_LMAX8_left.mat';
        ipsRoi     = 'LIPS0_cleaned.mat';
        mtRoi      = 'LTO1_O2_cleaned.mat';
        axisLims   = [-67 -18 -110 -40 -18 80];   
        wbAxLims   = [-86 0  -105 66 -2 60 ];

        vw         = [-75,30];
        slices     = {[-18 0 0],[0 -40 0],[0 0 -14 ]};
        lght       = 'left';
        fgUnculled = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep1/FP_LH_150_B2000_LMAX8.mat';

    case {'right'}
        feFileName = 'fe_culled_FP_150_B2000_LMAX8_right.mat';
        ipsRoi     = 'RIPS0_cleaned.mat';
        mtRoi      = 'RTO1_O2_cleaned.mat';
        axisLims   = [18 67 -110 -40 -18 80];
        wbAxLims   = [0 86 -105 66 -2 60 ];
        vw         = [75,30];
        slices     = {[18 0 0],[0 -40 0],[0 0 -14 ]};
        lght       = 'right';
        fgUnculled = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/life_mrtrix_rep1/FP_RH_150_B2000_LMAX8.mat';

    otherwise
        keyboard
end

% High-resolution Anatomy
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';
t1     = niftiRead(t1File);

% Load the FE structure
feDir      = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/fe_hemispheres';
fprintf('[%s] Loading the FE structure...\n',mfilename)
load(fullfile(feDir,feFileName));

% Load the ROIs
roiDir = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/mt_roi';
mtFileName   = fullfile(roiDir,mtRoi);
ips0FileName = fullfile(roiDir,ipsRoi);
mt   = dtiReadRoi(mtFileName);
ips0 = dtiReadRoi(ips0FileName);

% Find the fascicles in the connectome that touch both MT+ and IPS0.
operations = {'and','and'};
rois       = {mt,ips0};
fg         = feGet(fe,'fg acpc');
[mtIpsFG, keepFascicles] = feSegmentFascicleFromConnectome(fg, rois, operations, 'mt_ips_zero');
%[figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(mtIpsFG.fibers);

% Clean the fibers by length, fibers that too long are likely to go far
% frontal and not just touch MT+ and IPS0.
[Lnorm, Lmm] = mbaComputeFiberLengthDistribution(mtIpsFG, false);
maxSD        = 1; % Max standard deviation of the fibers to keep in the group.
fibers2delete  = Lnorm > maxSD;

% Now let's get the indices of the fibers in the FE structure:
fasIndices = find(keepFascicles);
fasIndices = fasIndices(fibers2delete);

% Now let's mark the fascicles as deleted.
keepFascicles(fasIndices) = false;

% Show the tract conencitng MT+ and IPS0
mtIpsFG = fgExtract(mtIpsFG,Lnorm < maxSD,'keep');
figName = sprintf('Test_MT_IPS0_connection_brain_%s_%s',hemisphere,feFileName(1:end-4));
figureHandle = mrvNewGraphWin(figName);
h  = mbaDisplayBrainSlice(t1, slices{1});
hold on
h  = mbaDisplayBrainSlice(t1, slices{2});
h  = mbaDisplayBrainSlice(t1, slices{3});
[figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(mtIpsFG.fibers,figureHandle,tractColor,'uniform');
delete(lightHandle)
view(vw(1),vw(2));
axis(axisLims);
lightHandle = camlight(lght);
set(gcf,'Position',[0.0148 0.0148 .35 .87])
drawnow
saveFig(figureHandle,fullfile(saveDir,figName),'jpg')

% Show all the fibers conecting MT before and after culling.
% After culling
operations = {'and'};
rois       = {mt};
fg         = feGet(fe,'fg acpc');
[mtAfterFG, keepFascicles] = feSegmentFascicleFromConnectome(fg, rois, operations, 'mt_after');

figName = sprintf('MT_connection_After_brain_%s_%s',hemisphere,feFileName(1:end-4));
figureHandle = mrvNewGraphWin(figName);
h  = mbaDisplayBrainSlice(t1, slices{1});
hold on
%h  = mbaDisplayBrainSlice(t1, slices{2});
h  = mbaDisplayBrainSlice(t1, slices{3});
[figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(mtAfterFG.fibers,figureHandle,tractColor,'uniform');
delete(lightHandle)
view(vw(1),vw(2));
axis(wbAxLims);
lightHandle = camlight(lght);
set(gcf,'Position',[0.0148 0.0148 .35 .87])
drawnow
saveFig(figureHandle,fullfile(saveDir,figName),'jpg')

% BEfore culling
fg = fgRead(fgUnculled);

% Find the fascicles in the connectome that touch MT+.
operations = {'and'};
rois       = {mt};
[mtBeforeFG, keepFascicles] = feSegmentFascicleFromConnectome(fg, rois, operations, 'mt_before');
figName = sprintf('MT_connection_Before_brain_%s_%s',hemisphere,feFileName(1:end-4));
figureHandle = mrvNewGraphWin(figName);
h  = mbaDisplayBrainSlice(t1, slices{1});
hold on
%h  = mbaDisplayBrainSlice(t1, slices{2});
h  = mbaDisplayBrainSlice(t1, slices{3});
[figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(mtBeforeFG.fibers,figureHandle,tractColor,'uniform');
delete(lightHandle)
view(vw(1),vw(2));
axis(wbAxLims);
lightHandle = camlight(lght);
set(gcf,'Position',[0.0148 0.0148 .35 .87])
drawnow
saveFig(figureHandle,fullfile(saveDir,figName),'jpg')

end % Main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName,type)

% MAke sure the folder to save the figure exists
[p,f,e] = fileparts(figName);
[success,message] = mkdir(p);
if ~isempty(message), disp(sprintf('%s.',message));end

% Find out which type of figure and geenerate the proper printing command.
switch type
    case 'eps'
        printCommand = sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName);
    case 'jpg'
        printCommand = sprintf('print(%s,  ''-djpeg'',''-r500'' , ''-noui'', ''%s'');', num2str(h),figName);
    otherwise
        keyboard
end

% do the printing here:
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
eval(printCommand);
end