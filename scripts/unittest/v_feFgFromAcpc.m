function v_feFgFromAcpc
% v_feFgFromAcpc.m
%
% Tests the feSet call that loads a fiber group from disk
% 
%

fe = feCreate;
dwiFile = fullfile(feDataPath('diffusion'),'life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz');
tempNi = niftiRead(dwiFile);
fe = feSet(fe, 'img2acpc xform', tempNi.qto_xyz);
fe = feSet(fe, 'acpc2img xform', inv(tempNi.qto_xyz));
fg = fullfile(feDataPath('tractography'),'life_demo_mrtrix_csd_lmax10_deterministic.mat');
fe = feSet(fe,'fg from acpc',fg);
assert(size(fe.fg.fibers,1)==50000,'Wrong number of fibers');
keyboard