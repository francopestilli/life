function s_ms_batch_multi_test_proclus(runType)

% Select the parameters for the current conditions
[algo, lmax, diffModParamsType] = getConditions(runType);

% Run the condition
s_ms_test_connectomes(algo,lmax,diffModParamsType);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [algo, lmax, diffModParamsType] = getConditions(runType)


switch runType
  %% MRTRIX probabilistic tractogrpahy

  % Beherens (1,0)
  case 1
    algo = 'p';lmax = 2;diffModParamsType = [1,0];
  case 2
    algo = 'p';lmax = 4;diffModParamsType = [1,0];
  case 3
    algo = 'p';lmax = 6;diffModParamsType = [1,0];
  case 4
    algo = 'p';lmax = 8;diffModParamsType = [1,0];
  case 5
    algo = 'p';lmax = 10;diffModParamsType = [1,0];
  case 6
    algo = 'p';lmax = 12;diffModParamsType = [1,0];
  
    % Optimal tensor form Corpus Callosum
  case 7
    algo = 'p';lmax = 2; diffModParamsType = [1.8,0.2];
  case 8
    algo = 'p';lmax = 4; diffModParamsType = [1.8,0.2];
  case 9
    algo = 'p';lmax = 6; diffModParamsType = [1.8,0.2];
  case 10
    algo = 'p';lmax = 8; diffModParamsType = [1.8,0.2];
  case 11
    algo = 'p';lmax = 10;diffModParamsType = [1.8,0.2];
  case 12
    algo = 'p';lmax = 12;diffModParamsType = [1.8,0.2];
   
   % isotropic tensor (1,0.9)
  case 1+12
    algo = 'p';lmax = 2;diffModParamsType = [1,0.9];
  case 2+12
    algo = 'p';lmax = 4;diffModParamsType = [1,0.9];
  case 3+12
    algo = 'p';lmax = 6;diffModParamsType = [1,0.9];
  case 4+12
    algo = 'p';lmax = 8;diffModParamsType = [1,0.9];
  case 5+12
    algo = 'p';lmax = 10;diffModParamsType = [1,0.9];
  case 6+12
    algo = 'p';lmax = 12;diffModParamsType = [1,0.9];
  
    % Small tensor
  case 7+12
    algo = 'p';lmax = 2; diffModParamsType = [1,0.1];
  case 8+12
    algo = 'p';lmax = 4; diffModParamsType = [1,0.1];
  case 9+12
    algo = 'p';lmax = 6; diffModParamsType = [1,0.1];
  case 10+12
    algo = 'p';lmax = 8; diffModParamsType = [1,0.1];
  case 11+12
    algo = 'p';lmax = 10;diffModParamsType = [1,0.1];
  case 12+12
    algo = 'p';lmax = 12;diffModParamsType = [1,0.1];
       
    % Large tensor
  case 1+24
    algo = 'p';lmax = 2; diffModParamsType = [1,0.3];
  case 2+24
    algo = 'p';lmax = 4; diffModParamsType = [1,0.3];
  case 3+24
    algo = 'p';lmax = 6; diffModParamsType = [1,0.3];
  case 4+24
    algo = 'p';lmax = 8; diffModParamsType = [1,0.3];
  case 5+24
    algo = 'p';lmax = 10;diffModParamsType = [1,0.3];
  case 6+24
    algo = 'p';lmax = 12;diffModParamsType = [1,0.3];
      
    % Small tensor canonical
  case 7+24
    algo = 'p';lmax = 2; diffModParamsType = [1,0.6];
  case 8+24
    algo = 'p';lmax = 4; diffModParamsType = [1,0.6];
  case 9+24
    algo = 'p';lmax = 6; diffModParamsType = [1,0.6];
  case 10+24
    algo = 'p';lmax = 8; diffModParamsType = [1,0.6];
  case 11+24
    algo = 'p';lmax = 10;diffModParamsType = [1,0.6];
  case 12+24
    algo = 'p';lmax = 12;diffModParamsType = [1,0.6];

    
  %% Deterministic tractography
  % Beherens (1,0)
  case 1+36
    algo='d';lmax = 2;diffModParamsType = [1,0];
  case 2+36
    algo='d';lmax = 4;diffModParamsType = [1,0];
  case 3+36
    algo='d';lmax = 6;diffModParamsType = [1,0];
  case 4+36
    algo='d';lmax = 8;diffModParamsType = [1,0];
  case 5+36
    algo='d';lmax = 10;diffModParamsType = [1,0];
  case 6+36
    algo='d';lmax = 12;diffModParamsType = [1,0];
 
    % Optimal tensor form Corpus Callosum
  case 7+36
    algo='d';lmax = 2; diffModParamsType = [1.8,0.2];
  case 8+36
    algo='d';lmax = 4; diffModParamsType = [1.8,0.2];
  case 9+36
    algo='d';lmax = 6; diffModParamsType = [1.8,0.2];
  case 10+36
    algo='d';lmax = 8; diffModParamsType = [1.8,0.2];
  case 11+36
    algo='d';lmax = 10;diffModParamsType = [1.8,0.2];
  case 12+36
    algo='d';lmax = 12;diffModParamsType = [1.8,0.2];
   
   % isotropic tensor (1,0.9)
  case 1+12+36
    algo='d';lmax = 2;diffModParamsType = [1,0.9];
  case 2+12+36
    algo='d';lmax = 4;diffModParamsType = [1,0.9];
  case 3+12+36
    algo='d';lmax = 6;diffModParamsType = [1,0.9];
  case 4+12+36
    algo='d';lmax = 8;diffModParamsType = [1,0.9];
  case 5+12+36
    algo='d';lmax = 10;diffModParamsType = [1,0.9];
  case 6+12+36
    algo='d';lmax = 12;diffModParamsType = [1,0.9];

    % Small tensor
  case 7+12+36
    algo='d';lmax = 2; diffModParamsType = [1,0.1];
  case 8+12+36
    algo='d';lmax = 4; diffModParamsType = [1,0.1];
  case 9+12+36
    algo='d';lmax = 6; diffModParamsType = [1,0.1];
  case 10+12+36
    algo='d';lmax = 8; diffModParamsType = [1,0.1];
  case 11+12+36
    algo='d';lmax = 10;diffModParamsType = [1,0.1];
  case 12+12+36
    algo='d';lmax = 12;diffModParamsType = [1,0.1];
       
    % Large tensor
  case 1+24+36
    algo='d';lmax = 2; diffModParamsType = [1,0.3];
  case 2+24+36
    algo='d';lmax = 4; diffModParamsType = [1,0.3];
  case 3+24+36
    algo='d';lmax = 6; diffModParamsType = [1,0.3];
  case 4+24+36
    algo='d';lmax = 8; diffModParamsType = [1,0.3];
  case 5+24+36
    algo='d';lmax = 10;diffModParamsType = [1,0.3];
  case 6+24+36
    algo='d';lmax = 12;diffModParamsType = [1,0.3];
       
    % Small tensor canonical with the same ration md/rd (1.8/0.2=1/.11)
  case 7+24+36
    algo='d';lmax = 2; diffModParamsType = [1,0.6];
  case 8+24+36
    algo='d';lmax = 4; diffModParamsType = [1,0.6];
  case 9+24+36
    algo='d';lmax = 6; diffModParamsType = [1,0.6];
  case 10+24+36
    algo='d';lmax = 8; diffModParamsType = [1,0.6];
  case 11+24+36
    algo='d';lmax = 10;diffModParamsType = [1,0.6];
  case 12+24+36
    algo='d';lmax = 12;diffModParamsType = [1,0.6];
    
  otherwise
    keyboard
end
  end