function wmMask = feMakeWMmask(classFile,mdFile,wmMaskFileName,prctThr)
% 
% Generates a White-matter mask using a freesurfer class file and the 
% mean diffusivity in diffusion data set. '
%
% The mask will contain voxels classigfed as white matter and will not
% contain gray matter nor CSF and ventricles. White-matter voxels are
% indicated by a 1 in the mask.
%  
%    wmMaskNifti = feMakeWMmask(classNifti,mdNifti,[wmMaskNiftiFname],[prctThr]);
% 
% INPUTS:
%     classFile - This is a 'classification' file. For example as generated
%                 using FreeSurfer. The is expected to be anifti or a path-to-a-nifti
%                 The nifti should be at the resolution fo the Anatomy volume 
%                 (e.g., 1mm or 0.7mm etc). The image values are expected to be:
%                 - The white matter is indicated as 3 and 4 (left and right
%                   respectively).
%                 - The gray  matter is classified as 5 and 6 (left and
%                   right respectively).
%        mdFile - This is either a path-to-a-nifti file or a dt6 file as
%                 generated in mrDiffusion. The nifti file should contain a
%                 3D volume with the mean diffusivity in the third dimension. 
%                 If a dt6 file is passed in the mean diffusivity image will
%                 be generated from the information contained in the dt6
%                 file.
% wmMaskFileName- This is the full path tot the nifti file that needs to be
%                 saved to disk. If this is not passed the fiel is not saved.
%        prctThr- Percentile (values 0-100) to use as threshold for the MD. 
%                 This parameter will define the lagest MD value accepted as 
%                 white matter, values higher than this percentile will be 
%                 classified as CSF and will not be part of the final WM mask.
%                 The default value (98) seem to work to identify the
%                 ventricles in the tested dataset.
%        
% OUTPUTS:
%    wmMask - This is a nifti structure containing the final white-matter
%             mask, resampled at the resolution fo the diffusion data.
%
% NOTES:
%    We  perform three operations. 
%      1. We create a mean diffusivity image from the tensor fit in a dt6 file. 
%      2. We create a white-matter mask forma classification file created
%         by freesurfer.
%      3. We restrict the white matter mask to only the voxels with mad
%         diffusivity below a certain value, white-matter has low MD values. This
%         value is take fromt he distribution fo the MD in the diffusion data set.
% 
%    The following is a related article to the method used here to generate a 
%    White-matter mask:
%      Jeurissen B, Leemans A, Tournier, J-D, Jones, DK and Sijbers J (2012).
%      Investigating the prevalence of complex fiber configurations in white
%      matter tissue with diffusion magnetic resonance imaging. Human Brain
%      Mapping doi: 10.1002/hbm.2209
%
% Franco (c) Stanford Vista Team 2012

% Check inputs
if (notDefined('classFile') || ~exist(classFile,'file'))
   error('[%s] The full path to a freesurfer class file is necessary.',mfilename);
else
  % Load a T1 class file created by freesurfer
  classFile = niftiRead(classFile);
end
if notDefined('prctThr'), prctThr = 98;end

% For the mean-diffusivity file we expect either the path to a dt6File or
% that to a nifti file with precomputed mean diffusivity
if ischar(mdFile)
  [~,~,e] = fileparts(mdFile);
  if strcmpi(e,'.mat')
    % A path to a dt6File was passed.
    % We compute a mean diffusivity nifti.
    tempFileName = sprintf('%s.nii.gz',tempname);
    [~, mdFile]  = dtiComputeMeanDiffusivity(mdFile,tempFileName);
  else
    % We assume that the path to a nfiti file of mean diffusivity was
    % passed in
    mdFile = niftiRead(mdFile);
  end

else
  error('[%s] Either the path to a dt6File or that to a MD nifti file is necessary...', mfilename) 
end

% Preprocess the class file from freesurfer to generate a white-matter
% mask.
%
% The white matter is classified as 3 and 4 (left/right) by freesurfer
% The gray  matter is classified as 5 and 6 (left/right)
%
% We set everything thatis NOT white-matter to 0.
classFile.data(classFile.data == 1) = 0;
classFile.data(classFile.data == 5) = 0;
classFile.data(classFile.data == 6) = 0; % CSF

% We set the white-matter to 1
classFile.data(classFile.data == 3) = 1;
classFile.data(classFile.data == 4) = 1;
% Check the the image: makeMontage3(classFile.data);

% Downsample the white-matter mask just created to the resolution of the
% mean-diffusivity image. The final WM mask will be generated at this reolution.
wmMask       = mrAnatResampleToNifti(classFile, mdFile);

% Now let's remove from the mask the ventricles. To do so we remove voxels
% that have high mean diffusivity.
% 
% Find the indices of the MD imag that are above a certain trheshold
notcsf = false(size(mdFile.data));
maxMD  = prctile(mdFile.data(:),prctThr); % We set a cutoff at the 96.5% percentile 
notcsf( (mdFile.data < maxMD) ) = 1;

% Now only accept the voxels that were in the original white-matter
% segmentation but that are not classified as CSF given their mean
% diffusivity value.
wmMask.data = single(logical(wmMask.data) & notcsf);

% Save the new white-matter mask to disk:
if ~notDefined('wmMaskFileName')
  fprintf('[%s] Saving WM mask to disk: %s..\n.',mfilename,wmMaskFileName)
  wmMask.fname   = wmMaskFileName;
  wmMask.descrip = [which(mfilename),'- Original File: ',classFile.fname];
  niftiWrite(wmMask);
end

return

