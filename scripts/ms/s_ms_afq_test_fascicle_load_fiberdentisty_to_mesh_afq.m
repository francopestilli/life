% Example of how to make a fider density map fr om a faccile, and render it
% on a a surface.

clx

%cd  /azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/fe_arcuate_cst_test
cd /azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/fe_arcuate_importance/
% Prepare the path for the nifti to the0 segemntation file
cortex = '/biac2/wandell2/data/anatomy/pestilli/t1_class_twovalued.nii.gz';
thresh = 0.0000125; % Threshold for the overlay image
crange = [.0000125 .5]; % Color range of the overlay image


fiberDensityName = 'SLF_left_fiberDensity_tensor_best_fibers';

% Load the fe structure:
%load arcuateCstUnionFE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax12_m_prob-500000_diffModAx100Rd0_arcuateCstUnionFE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax12_m_prob-500000_diffModAx100Rd0_.mat
load  FE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_dwi_0005_0_tensor-500000_diffModAx100Rd0_lSlf_sd300_lSlfFE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_dwi_0005_0_tensor-500000_diffModAx100Rd0_lSlf_sd300_lSlf.mat
%load  FE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_dwi_0005_0_tensor-500000_diffModAx100Rd0_rSlf_sd300_rSlfFE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_dwi_0005_0_tensor-500000_diffModAx100Rd0_rSlf_sd300_rSlf.mat

% Make a map of the fiber density form the fascicle.
feMakeFiberDensityNifti(fe, fiberDensityName, 11,1);

% Render the cortical surface colored by the arcuate endpoint density 
[p msh]= AFQ_RenderCorticalSurface(cortex, 'overlay' , fiberDensityName, 'crange', crange, 'thresh', thresh,'cmap','autumn');

% Work on the figure
set(gcf,'color','k')
axis off
view(-90,0); %camlight
eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(gcf),sprintf('~/Dropbox/%s',fiberDensityName)));

% Make a figure of the fascicle on a brain slice
fg = feGet(fe,'fibers acpc');
w  = feGet(fe,'fiber weights');
fg.fibers = fg.fibers(w > 0.002);

feConnectomeDisplay( feSplitLoopFibers( fg ),figure,[.6 .23 .4],[],[],.1);
view(-90,0); 
hold on
t1 = '/biac2/wandell2/data/anatomy/pestilli/t1.nii.gz';
mctDisplayBrainSlice(niftiRead( t1 ),[-10 0 0])
camlight
eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(gcf),sprintf('~/Dropbox/%s_fascicle_on_brain',fiberDensityName)));


fiberDensityName = 'SLF_left_fiberDensity_probabilistic_best_fibers';

% Load the fe structure:
%load arcuateCstUnionFE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax2__stream-500000_diffModAx100Rd0_arcuateCstUnionFE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax2__stream-500000_diffModAx100Rd0_.mat
%load FE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax8__m_prob-500000_diffModAx100Rd0_lSlf_sd300_lSlfFE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax8__m_prob-500000_diffModAx100Rd0_lSlf_sd300_lSlf.mat
load FE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax8__m_prob-500000_diffModAx100Rd0_rSlf_sd300_rSlfFE_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax8__m_prob-500000_diffModAx100Rd0_rSlf_sd300_rSlf.mat


% Make a map of the fiber density form the fascicle.
feMakeFiberDensityNifti(fe, fiberDensityName, 11,1);

% Render the cortical surface colored by the arcuate endpoint density 
p = AFQ_RenderCorticalSurface(cortex, 'overlay' , fiberDensityName, 'crange', crange, 'thresh', thresh,'cmap','autumn');

% Work on the figure
set(gcf,'color','k')
axis off
view(-90,0); %camlight
%set(gcf,'Position',[0 0 .45 .95],'Color',[0 0 0]);  drawnow
eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(gcf),sprintf('~/Dropbox/%s',fiberDensityName)));

% Make a figure of the fascicle on a brain slice
fg = feGet(fe,'fibers acpc');
w  = feGet(fe,'fiber weights');
fg.fibers = fg.fibers(w > 0.002);

feConnectomeDisplay( feSplitLoopFibers( fg ),figure,[.6 .23 .4],[],[],.1);
view(-90,0); 
hold on
t1 = '/biac2/wandell2/data/anatomy/pestilli/t1.nii.gz';
mctDisplayBrainSlice(niftiRead( t1 ),[-10 0 0])
camlight
eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(gcf),sprintf('~/Dropbox/%s_fascicle_on_brain',fiberDensityName)));

% % make a histogram with fibers before and after life.
%fdNiiall  = feMakeFiberDensityNifti(fe, fiberDensityName, 0,0);
%fdPreLife = fdNiiall.data(fdNiiall.data(:) ~= 0);

% % Make a map of the fiber density form the fascicle.
%fdNiiBest = feMakeFiberDensityNifti(fe, fiberDensityName, 0,1);
%fdPostLife = fdNiiBest.data(fdNiiBest.data(:) ~= 0);
