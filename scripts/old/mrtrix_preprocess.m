% () in matlab
baseDir    = fullfile(mrvDataRootPath,'diffusion','sampleData');
dtFile     = fullfile(baseDir,'dti40','dt6.mat');
dwiFile    = fullfile(baseDir,'raw','dwi.nii.gz');

% (1) in terminal:
cd(fullfile(baseDir,'dti40')) 

% (2) in matlab
mkdir mrtrix
% (3) in matlab
cd mrtrix;
gunzip -c ../bin/brainMask.nii.gz > mask.nii

% (4) in shell
mrconvert mask.nii mask.mif

% (5) in MATLAB:
grads=[dlmread('../../raw2/dti30_aligned_trilin.bvecs'); dlmread('../../raw2/dti30_aligned_trilin.bvals')];
% check bvecs & bvals
bvecs=dlmread('../../raw2/dti30_aligned_trilin.bvecs'); bvals = dlmread('../../raw2/dti30_aligned_trilin.bvals')
size(bvecs) % should be 3 70
size(bvals) % should be 1 70
grads = [bvecs;bvals] % create grads
size(grads) % should be 4 70
% (6) write grads
dlmwrite('grads',grads',' ')
% (7) in terminal
gunzip -c  ../../raw2/dti30_aligned_trilin.nii.gz > dti30_aligned_trilin.nii
% (8)
dwi2tensor -grad grads dti30_aligned_trilin.nii dt.mif
% (9) Compute a fa.mif file 
tensor2FA dt.mif - | mrmult - mask.mif fa.mif
% (10) Compute a pdd file 
tensor2vector dt.mif - | mrmult - fa.mif pdd.mif
% (11) View the pdd file
mrview pdd.mif
% (12) Create single-fiber mask 
erode mask.mif - | erode - - | mrmult fa.mif - - | threshold - -abs 0.7 sf.mif
% (13) Create response function coefficients (~30 sec) 
estimate_response -grad grads dti30_aligned_trilin.nii sf.mif response.txt
% (14) Constrained Spherical Deconvolution (~30 min)
csdeconv -grad grads dti30_aligned_trilin.nii response.txt -mask mask.mif CSD.mif
% or
csdeconv -grad grads dti30_aligned_trilin.nii response.txt -lmax 10 -mask mask.mif CSD10.mif
% (15) Generate a white matter mask for whole-brain tractography 
erode mask.mif - | erode - - | mrmult fa.mif - - | threshold - -abs 0.25 wm.mif
% (16) Probabilistic
streamtrack SD_PROB CSD10.mif -seed wm.mif -mask mask.mif all_100K.tck -num 100000 -trials 500
% or
streamtrack SD_PROB CSD10.mif -seed wm.mif -mask mask.mif all_1000K.tck -num 1000000 -trials 1000
% (17) in MATLAB.
fg = dtiImportFibersMrtrix('all_100K.tck');
% (18) Create pdb
mtrExportFibers(fg, 'all_100K.pdb', eye(4));