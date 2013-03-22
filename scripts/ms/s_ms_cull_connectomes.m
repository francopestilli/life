function feFileToSave = s_ms_cull_connectomes(connectomeType,rep)
%
% feFileToSave = s_ms_cull_connectomes(connectomeType,rep)
%
% Loads a connectome and culls it down to the minimum number of fibers.
%
% Written by Franco Pestilli (c) Stanford University 2013 
if notDefined('rep'), rep = [1,2,3];end

% We have three reps of each conenctome.
for irep = 1:length(rep)
    [trackingType, lmax, bval, diffusionModelParams] = getConditions(connectomeType);
    
    % Information on the path to the files to load.
    % This is where the inputs will be loaded from
    feFileToLoad = msBuildFeFileName(trackingType,lmax,bval,rep(irep),diffusionModelParams);
    
    % Get the fe structure
    fprintf('[%s] loading the LiFE structure...\n',mfilename);
    if exist(feFileToLoad,'file')
        fprintf('[%s] File exist OK.\n',mfilename)
    else
        keyboard
    end
    load(feFileToLoad);
    
    fprintf('[%s] Culling LiFE structure...\n',mfilename);
    [fe, cullingInfo] = feConnectomeCull(fe,1000,'sgdnn');
    
    feFileToSave = [feFileToLoad(1:end-4),'culledL2','.mat'];
    fprintf('[%s] Saving the culled LiFE structure...\n%s\n',mfilename,feFileToSave);
    feFileToSave = [feFileToLoad(1:end-4),'culledL2','.mat'];
    save(feFileToSave,'fe','cullingInfo','-v7.3');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [algo, lmax, bval, diffModParamsType] = getConditions(runType)
%
% trackingType,lmax,bval,rep,diffusionModelParams
%
% Returns the desired options for a tractography.
switch runType
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 2000
    case 1
        algo = 'p';lmax = 2;bval=2000;diffModParamsType = [1,0];
    case 2
        algo = 'p';lmax = 4;bval=2000;diffModParamsType = [1,0];
    case 3
        algo = 'p';lmax = 6;bval=2000;diffModParamsType = [1,0];
    case 4
        algo = 'p';lmax = 8;bval=2000;diffModParamsType = [1,0];
    case 5
        algo = 'p';lmax = 10;bval=2000;diffModParamsType = [1,0];
    case 6
        algo = 'p';lmax = 12;bval=2000;diffModParamsType = [1,0];
    case 7
        algo = 'p';lmax = 14;bval=2000;diffModParamsType = [1,0];
    case 8
        algo = 'p';lmax = 16;bval=2000;diffModParamsType = [1,0];
      
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 4000
    case 9
        algo = 'p';lmax = 2;bval=4000;diffModParamsType = [1,0];
    case 10
        algo = 'p';lmax = 4;bval=4000;diffModParamsType = [1,0];
    case 11
        algo = 'p';lmax = 6;bval=4000;diffModParamsType = [1,0];
    case 12
        algo = 'p';lmax = 8;bval=4000;diffModParamsType = [1,0];
    case 13
        algo = 'p';lmax = 10;bval=4000;diffModParamsType = [1,0];
    case 14
        algo = 'p';lmax = 12;bval=4000;diffModParamsType = [1,0];
    case 15
        algo = 'p';lmax = 14;bval=4000;diffModParamsType = [1,0];
    case 16
        algo = 'p';lmax = 16;bval=4000;diffModParamsType = [1,0];
         
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 1000
    case 17
        algo = 'p';lmax = 2;bval=1000;diffModParamsType = [1,0];
    case 18
        algo = 'p';lmax = 4;bval=1000;diffModParamsType = [1,0];
    case 19
        algo = 'p';lmax = 6;bval=1000;diffModParamsType = [1,0];
    case 20
        algo = 'p';lmax = 8;bval=1000;diffModParamsType = [1,0];
    case 21
        algo = 'p';lmax = 10;bval=1000;diffModParamsType = [1,0];
    case 22
        algo = 'p';lmax = 12;bval=1000;diffModParamsType = [1,0];
    case 23
        algo = 'p';lmax = 14;bval=1000;diffModParamsType = [1,0];
    case 24
        algo = 'p';lmax = 16;bval=1000;diffModParamsType = [1,0];
         
    % MRTRIX deterministic tractogrpahy 
    % BVAL = 2000
    case 1+16
        algo = 'd';lmax = 2;bval=2000;diffModParamsType = [1,0];
    case 2+16
        algo = 'd';lmax = 4;bval=2000;diffModParamsType = [1,0];
    case 3+16
        algo = 'd';lmax = 6;bval=2000;diffModParamsType = [1,0];
    case 4+16
        algo = 'd';lmax = 8;bval=2000;diffModParamsType = [1,0];
    case 5+16
        algo = 'd';lmax = 10;bval=2000;diffModParamsType = [1,0];
    case 6+16
        algo = 'd';lmax = 12;bval=2000;diffModParamsType = [1,0];
    case 7+16
        algo = 'd';lmax = 14;bval=2000;diffModParamsType = [1,0];
    case 8+16
        algo = 'd';lmax = 16;bval=2000;diffModParamsType = [1,0];
      
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 4000
    case 9+16
        algo = 'd';lmax = 2;bval=4000;diffModParamsType = [1,0];
    case 10+16
        algo = 'd';lmax = 4;bval=4000;diffModParamsType = [1,0];
    case 11+16
        algo = 'd';lmax = 6;bval=4000;diffModParamsType = [1,0];
    case 12+16
        algo = 'd';lmax = 8;bval=4000;diffModParamsType = [1,0];
    case 13+16
        algo = 'd';lmax = 10;bval=4000;diffModParamsType = [1,0];
    case 14+16
        algo = 'd';lmax = 12;bval=4000;diffModParamsType = [1,0];
    case 15+16
        algo = 'd';lmax = 14;bval=4000;diffModParamsType = [1,0];
    case 16+16
        algo = 'd';lmax = 16;bval=4000;diffModParamsType = [1,0];
         
    % MRTRIX probabilistic tractogrpahy 
    % BVAL = 1000
    case 17+16
        algo = 'd';lmax = 2;bval=1000;diffModParamsType = [1,0];
    case 18+16
        algo = 'd';lmax = 4;bval=1000;diffModParamsType = [1,0];
    case 19+16
        algo = 'd';lmax = 6;bval=1000;diffModParamsType = [1,0];
    case 20+16
        algo = 'd';lmax = 8;bval=1000;diffModParamsType = [1,0];
    case 21+16
        algo = 'd';lmax = 10;bval=1000;diffModParamsType = [1,0];
    case 22+16
        algo = 'd';lmax = 12;bval=1000;diffModParamsType = [1,0];
    case 23+16
        algo = 'd';lmax = 14;bval=1000;diffModParamsType = [1,0];
    case 24+16
        algo = 'd';lmax = 16;bval=1000;diffModParamsType = [1,0];
        
    % MRTRIX tensor-based tractogrpahy 
    % BVAL = 2000,4000,1000
    case 25+16
        algo = 't';lmax = 2;bval=2000;diffModParamsType = [1,0];
    case 26+16
        algo = 't';lmax = 2;bval=4000;diffModParamsType = [1,0];
    case 27+16
        algo = 't';lmax = 2;bval=1000;diffModParamsType = [1,0];
  
  otherwise
        keyboard
end

end
