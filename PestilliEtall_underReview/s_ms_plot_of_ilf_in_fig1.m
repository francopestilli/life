% Move to the folder where the connectome structures are saved
cd /azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/fe_structures/

% THis is how I created the conenctomes
% Load a high-lmax probabilistic fe to get wigly fibers
%load 0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax2_0_stream-500000_diffModAx100Rd0_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax2_0_stream-500000_diffModAx100Rd0_.mat
%fg = feSplitLoopFibers( feGet(fe,'fg acpc') );
%fgWrite(fg,'fg_streamLmax2Split','mat');

% Load a high-lmax probabilistic fe to get wigly fibers
%load 0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax10__m_prob-500000_diffModAx100Rd0_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax10__m_prob-500000_diffModAx100Rd0_.mat
%fg = feSplitLoopFibers( feGet(fe,'fg acpc') );
%fgWrite(fg,'fg_probLmax10Split','mat');

% Load the ILF segments created from the above connectomes.
cd Example_tracts_from_occipital_connectomes/fascicles/

fg = fgRead('mrtrix_lmax12_prob_ILF_right_hemisphere.mat');

% Remove the fibers outliers
%fg =AFQ_removeFiberOutliers(fg,3,40,100);

% Display the fascicle
feConnectomeDisplay( feSplitLoopFibers(fg), figure, [.85 .6 .7], [], [], .25);
view(30,0)
eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', ...
  num2str(gcf),'mrtrix_lmax12_prob_ILF_right_hemisphere_view30_0.jpg'));

% Load the Lmax=2 fascicle
fg = fgRead('mrtrix_lmax2_stream_ILF_rigth_hemisphere.mat');

% Remove the fibers outliers
%fg =AFQ_removeFiberOutliers(fg,3.5,40,100);

% Display the fascicle
feConnectomeDisplay( feSplitLoopFibers(fg), figure, [.85 .6 .7], [], [], .25);
view(30,0)
eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', ...
  num2str(gcf),'mrtrix_lmax2_stream_ILF_rigth_hemisphere_view30_0.jpg'));

