function p = msPaths(pathType)
%
% Returns the paths of interest for the LiFE manuscript.
%
% Example: p = msCodePaths('code')
%      
% 
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.

if notDefined('pathType'), error('[%s] pathTYpe is requried as an input.',mfilename);end

switch pathType
  case 'projectDir'
    p = fullfile(msPaths('data'),'150dirs_b1000_b2000_b4000');
    
  case 'dt61000'
    p = fullfile(msPaths('session1000'),'150dirs_b1000_1');
    
  case 'dt62000'
    p = fullfile(msPaths('session2000'),'150dirs_b2000_1');
    
  case 'dt64000'
    p = fullfile(msPaths('session4000'),'150dirs_b4000_1');
    
  case 'session1000'
    p = fullfile(msPaths('projectDir'),'150dirs_b1000_b4000');
    
  case 'session2000'
    p = fullfile(msPaths('projectDir'),'150dirs_b2000');
    
  case 'session4000'
    p = fullfile(msPaths('s'),'150dirs_b1000_b4000');

  case {'code'}
    p = fileparts(which(mfilename));
    
  case {'data'}
    p = fullfile('/azure/scr1/frk');
    
  case {'bwrois','wernickerois','brocarois'}
    p = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/roi_broca_wernicke';
  case {'bwroisjw','wernickeroiswinawer','brocaroisjw'}
    p = '/azure/scr1/frk/JW_96dirs_b2000_1p5iso/results/life_mrtrix_rep1/broca_wernicke_roi';
  case {'mtrois'}
    p = '/azure/scr1/frk/JW_96dirs_b2000_1p5iso/results/life_mrtrix_rep1/mt_roi';
   
  otherwise
    error('[%s] Cannot find requested pathType (%s).',mfilename,pathType);
end