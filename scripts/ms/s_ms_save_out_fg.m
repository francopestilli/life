cd /azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/fe_structures

load 0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax2__m_prob-500000_diffModAx50Rd0_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax2__m_prob-500000_diffModAx50Rd0_.mat

fg = feSplitLoopFibers( feGet(fe,'fg img') );

fg = dtiXformFiberCoords(fg,feGet(fe,'xform 2 acpc'))

fg = rmfield(fg, 'pathwayInfo')

fgWrite(fg,['new2_FG'],'pdb')

feConnectomeDisplay(fg, figure)