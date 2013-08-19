function [feFileToLoad, feLoadName, feLoadDir] = msBuildFeFileName(trackingType,lmax,bval,rep,diffusionModelParams,cullType)
%
% Build a file name for one of the several connectomes preprocessed for the
% LiFE manuscript.
%  
%  [feFileToLoad, feLoadName] = msBuildFeFileName(trackingType,lmax,bval,rep,diffusionModelParams,cullType)
% 
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

% Build a file name for the directory to load the data
loadDir = fullfile( msPaths('projectDir'), 'results' );
if notDefined('cullType'), cullType='';end

% Build a file name for the FE file.
if ~any(lmax == [2,4,6,8,10,12,14,16])
  error('[%s] Lmax should one of the following: [2,4,6,8,10,12,14,16]\n',mfilename)
end

% bval files names have slightly different. So I organized them in order of bval.
switch bval
  case {1,1000,'1000'}
    bval = 1;
  case {2,2000,'2000'}
    bval = 2;
  case {4,4000,'4000'}
    bval = 3;
  otherwise
    keyboard
end

% Find the repetition number of trackrography
switch rep
  case {'one',1,'1'}
    connectSubfolders = 'life_mrtrix_rep1';
  case {'two',2,'2'}
    connectSubfolders = 'life_mrtrix_rep2';
  case {'three',3,'3'}
    connectSubfolders = 'life_mrtrix_rep3';
  otherwise
    keyboard
end

% Get the name of the connctome using lmax and the type of tractography
switch lower(trackingType)
  case {'deterministic','d','det'}
    connectomeFile = { ...
      sprintf( '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax%i_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_stream-500000.pdb',   lmax), ...
      sprintf( '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax%i_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_stream-500000.pdb',lmax),...
      sprintf( '0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_csd_lmax%i_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_brainmask_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_wm_stream-500000.pdb',   lmax),...
      };
  case {'probabilistic','p','prob'}
    connectomeFile = { ...
      sprintf( '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax%i_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_prob-500000.pdb',   lmax), ...
      sprintf( '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax%i_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb',lmax),...
      sprintf( '0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_csd_lmax%i_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_brainmask_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_wm_prob-500000.pdb',   lmax);
      };
    case {'tensor','t'}
    connectomeFile = { ...
          sprintf( '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_dwi_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_tensor-500000.pdb'), ...
          sprintf( '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_dwi_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_tensor-500000.pdb'),...
          sprintf( '0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_dwi_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_brainmask_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_wm_tensor-500000.pdb');
          };
otherwise
    keyboard
end

% Set up file names
% Find an identificative name for the connectome that is short enough:
cName        = [connectomeFile{bval}(1:57),'_',connectomeFile{bval}(end-16:end-4)];
feLoadDir    = fullfile(loadDir,connectSubfolders,'fe_structures');
feLoadName   = sprintf('%s_diffModAx%sRd%s_%s',cName,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)));
feFileToLoad = fullfile(feLoadDir,[feLoadName,feLoadName,cullType,'.mat']);

end