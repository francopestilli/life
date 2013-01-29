function p = msPaths(pathType)
%
% Returns the paths of interest for the LiFE manuscript.
%
% Example: p = msCodePaths('code')
%      
% 
% Franco (c) VISTA Team 2012

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
  
  otherwise
    error('[%s] Cannot find requested pathType (%s).',mfilename,pathType);
end