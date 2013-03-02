function s_ms_test_connectomes_maps(trackingType,dataType,lmax,diffusionModelParams,recompute)
%
% Load FE structeres obtained by preprocessing connectomesconstrained within a
% region of interest and within the cortex and makes some basic plot of
% statistics in the connectomes
%
% See also:
%
% Franco (C) 2012 Stanford VISTA team.

if notDefined('trackingType'),trackingType = 'tensor';end
if notDefined('lmax'),        lmax=[2 :2: 16];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('recompute'), recompute=1;end
if notDefined('dataType'), dataType='96dirs';end

% ROIs connectomes and saved paths
connectSubfolders = {'life_mrtrix_rep1'};
switch dataType
    case {'96dirs'}
        lmax = lmax(lmax<=12);
        bval = 2000;
        
    case {'150dirs'}
        lmax = lmax(lmax<=16);
        bval   =  [1000,2000,4000];
        
    otherwise
        keyboard
end

switch trackingType
  case {'t','tensor'}
    lmax=0;
end

for i_lmax = 1:length(lmax)
  switch dataType
    case {'150dirs','2mm'}
      
      % ROIs connectomes and saved paths
      projectDir           = '/azure/scr1/frk/150dirs_b1000_b2000_b4000';
      plots_saveDir     = fullfile('/home/frk/Dropbox','connectomes_mapshists_150dirs');
      if ~isdir(plots_saveDir), mkdir(plots_saveDir);end
      
      % Get the connectome names
      switch trackingType
        case {'deterministic','d'}
          connectomeFile = { ...
            sprintf( '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax%i_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_stream-500000.pdb',lmax(i_lmax)), ...
            sprintf( '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax%i_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_stream-500000.pdb',lmax(i_lmax)),...
            sprintf( '0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_csd_lmax%i_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_brainmask_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_wm_stream-500000.pdb',lmax(i_lmax)),...
            };
        case {'probabilistic','p'}
          connectomeFile = { ...
            sprintf( '0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_csd_lmax%i_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_brainmask_0009_01_DWI_2mm150dir_2x_b1000_aligned_trilin_wm_prob-500000.pdb',lmax(i_lmax)), ...
            sprintf( '0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_csd_lmax%i_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_brainmask_0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin_wm_prob-500000.pdb',lmax(i_lmax)),...
            sprintf( '0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_csd_lmax%i_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_brainmask_0005_01_DWI_2mm150dir_2x_b4000_aligned_trilin_wm_prob-500000.pdb',lmax(i_lmax));
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
    case {'96dirs','1.5mm'}
      
      % ROIs connectomes and saved paths
      projectDir           = '/azure/scr1/frk/FP_96dirs_b2000_1p5iso';
      plots_saveDir        = fullfile('/home/frk/Dropbox','connectomes_mapshists_fp96dirs');
      if ~isdir(plots_saveDir), mkdir(plots_saveDir);end
      
      % Get the connectome names
      switch trackingType
        case {'deterministic','d'}
          connectomeFile = { ...
            sprintf( 'run01_fliprot_aligned_trilin_csd_lmax%i_run01_fliprot_aligned_trilin_brainmask_run01_fliprot_aligned_trilin_wm_stream-500000.pdb',lmax(i_lmax))};
        case {'probabilistic','p'}
          connectomeFile = { ...
            sprintf( 'run01_fliprot_aligned_trilin_csd_lmax%i_run01_fliprot_aligned_trilin_brainmask_run01_fliprot_aligned_trilin_wm_prob-500000.pdb',lmax(i_lmax))};
        case {'tensor','t'}
          connectomeFile = { ...
            sprintf( 'run01_fliprot_aligned_trilin_dwi_run01_fliprot_aligned_trilin_brainmask_run01_fliprot_aligned_trilin_wm_tensor-500000.pdb')};
        otherwise
          keyboard
      end
    otherwise
      keyboard
  end
  
  % Location to  load the fe structure from
  loadDir  = fullfile(projectDir,'results');
  lmaxs    = num2str(lmax); lmaxs(isspace(lmaxs))='_';
  fileName = fullfile(plots_saveDir,sprintf('%s_ax%srd%s_tck%s_res_March2013_lmax%s.mat',mfilename,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)),upper(trackingType),lmaxs));
  
  if recompute % Recompute the satistics  
    % Loop over Bvalues and tracking reps.
    for i_bval = 1:length(connectomeFile)
      for irep = 1:length(connectSubfolders)
        
        % Set up file names
        % Find an identificative name for the connectome that is short enough:
        cName = [connectomeFile{i_bval}(1:57),'_',connectomeFile{i_bval}(end-16:end-4)];
        feLoadDir         = fullfile(loadDir,connectSubfolders{irep},'fe_structures');
        feLoadName        = sprintf('%s_diffModAx%sRd%s_%s',cName,num2str(100*diffusionModelParams(1)),num2str(100*diffusionModelParams(2)));
        
        % File to load
        feFileToLoad =  fullfile(feLoadDir,[feLoadName,feLoadName,'.mat']);
        
        if (exist(feFileToLoad,'file') == 2)
          fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)
          load(feFileToLoad);
        else
          fprintf('[%s] FE file NOT found: \n%s\n ======================================== \n\n\n',mfilename,feFileToLoad)
          keyboard
        end
        
        % Make a plot of the maps
        slice = 26:35;
        for is = 1:length(slice)
        [~, figH] = fePlot(fe,'fiberdensitymap',slice(is));
        for fi = 1:length(figH)
            figN = get(figH(fi),'name');
            figName = sprintf('%s_slc%i_lmx%s',figN(~isspace(figN)),slice(is),num2str(lmaxs(i_lmax)));
            
            % Save the figure
            saveFig(figH(fi),fullfile(plots_saveDir,figName))
        end
        end
        
        % Make a plot of the hists
        [~, figH] = fePlot(fe,'fiberdensityhist');
        for fi = 1:length(figH)
            figN = get(figH(fi),'name');
            figName = sprintf('%s_lmx%s',figN(~isspace(figN)),num2str(lmaxs(i_lmax)));
            
            % Save the figure
            saveFig(figH(fi),fullfile(plots_saveDir,figName))
        end
        close all; drawnow
      end % irep
    end   % i_bval
  end
end     % i_lmax

end % Main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)

printCommand = sprintf('print(%s, ''-djpeg90'', ''-opengl'', ''%s'')', num2str(h),figName);
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

% do the printing here:
eval(printCommand);

% save a wiki-compatible file
eval(sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName));

end