function s_ms_afq_arcuate_display_JW
%
% Uses AFQ to segment a connectome and generate 20 major faascicles. 
%
% Then extract the left and right arcuate fasciculum and display them on a slice.
%
% Franco (c) Stanford Vista Team 2012
  

% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

fas_colors = [.88 .78 .45];

algo = {'t','p'};
fascicleDir = '/azure/scr1/frk/JW_96dirs_b2000_1p5iso/results/life_mrtrix_rep1/fas_right_arcuate';
fascicleIndices = [19,20]; % Left and right arcuate

anatomy = '/biac2/wandell2/data/diffusion/winawer/20120410_2202/t1/t1.nii.gz';

for ia = 1:length(algo)
switch algo{ia}
    case {'t','tensor'}
        fasName = 'fas_run01_fliprot_aligned_trilin_dwi_run01_fliprot_aligned_tr_tensor-500000_diffModAx100Rd0_lArc_lArc_cleaned_sd300_lArc.mat';
    case {'p','prob'}
        fasName = 'fas_run01_fliprot_aligned_trilin_csd_lmax8_run01_fliprot_alig_m_prob-500000_diffModAx100Rd0_lArc_lArc_cleaned_sd300_lArc.mat'; 
    otherwise
        keyboard
end

fasPath     = fullfile(fascicleDir,fasName);
fprintf('[%s] Loading it: \n%s\n ======================================== \n\n',mfilename,fasPath)
load(fasPath);

% Show the fascicles
leftFG  = fascicles(fascicleIndices(1));
rightFG = fascicles(fascicleIndices(2));

% SAGITAL Right arcuate
figureHandle = figure;
feDisplayBrainSlice(niftiRead(anatomy), [28 0 0]);
hold on
[figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(rightFG.fibers,figureHandle,fas_colors,'single');
delete(lightHandle);
axis([-75 75 -75 57 -32 70]);
view(90,0);
lightHandle = camlight('right');
drawnow

% SAGITAL Left Arcuate
figureHandle = figure;
feDisplayBrainSlice(niftiRead(anatomy), [-28 0 0]);
hold on
[figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(leftFG.fibers,figureHandle,fas_colors,'single');
delete(lightHandle);
axis([-75 75 -75 57 -32 70]);
view(-90,0);
lightHandle = camlight('right');
drawnow;

% AXIAL Right arcuate
figureHandle = figure;
feDisplayBrainSlice(niftiRead(anatomy), [0 0 -4]);
hold on
[figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(rightFG.fibers,figureHandle,fas_colors,'single');
delete(lightHandle);

% AXIAL Left Arcuate
[figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(leftFG.fibers,figureHandle,fas_colors,'single');
delete(lightHandle);
axis([-75 75 -75 57 -32 70]);
view(0,90);
lightHandle = camlight('right');
drawnow;

end

end % End function

