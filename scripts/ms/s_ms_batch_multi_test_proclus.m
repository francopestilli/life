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
  case 7
    algo = 'p';lmax = 14;diffModParamsType = [1,0];
  case 8
    algo = 'p';lmax = 16;diffModParamsType = [1,0];
 
    % Optimal tensor form Corpus Callosum
  case 9
    algo = 'p';lmax = 2; diffModParamsType = [1.8,0.2];
  case 10
    algo = 'p';lmax = 4; diffModParamsType = [1.8,0.2];
  case 11
    algo = 'p';lmax = 6; diffModParamsType = [1.8,0.2];
  case 12
    algo = 'p';lmax = 8; diffModParamsType = [1.8,0.2];
  case 13
    algo = 'p';lmax = 10;diffModParamsType = [1.8,0.2];
  case 14
    algo = 'p';lmax = 12;diffModParamsType = [1.8,0.2];
  case 15
    algo = 'p';lmax = 14;diffModParamsType = [1.8,0.2];
  case 16
    algo = 'p';lmax = 16;diffModParamsType = [1.8,0.2];
   
   % isotropic tensor (1,1)
  case 1+16
    algo = 'p';lmax = 2;diffModParamsType = [1,1];
  case 2+16
    algo = 'p';lmax = 4;diffModParamsType = [1,1];
  case 3+16
    algo = 'p';lmax = 6;diffModParamsType = [1,1];
  case 4+16
    algo = 'p';lmax = 8;diffModParamsType = [1,1];
  case 5+16
    algo = 'p';lmax = 10;diffModParamsType = [1,1];
  case 6+16
    algo = 'p';lmax = 12;diffModParamsType = [1,1];
  case 7+16
    algo = 'p';lmax = 14;diffModParamsType = [1,1];
  case 8+16
    algo = 'p';lmax = 16;diffModParamsType = [1,1];
 
    % Small tensor
  case 9+16
    algo = 'p';lmax = 2; diffModParamsType = [.5,0];
  case 10+16
    algo = 'p';lmax = 4; diffModParamsType = [.5,0];
  case 11+16
    algo = 'p';lmax = 6; diffModParamsType = [.5,0];
  case 12+16
    algo = 'p';lmax = 8; diffModParamsType = [.5,0];
  case 13+16
    algo = 'p';lmax = 10;diffModParamsType = [.5,0];
  case 14+16
    algo = 'p';lmax = 12;diffModParamsType = [.5,0];
  case 15+16
    algo = 'p';lmax = 14;diffModParamsType = [.5,0];
  case 16+16
    algo = 'p';lmax = 16;diffModParamsType = [.5,0];
      
    % Large tensor
  case 1+32
    algo = 'p';lmax = 2; diffModParamsType = [2,0];
  case 2+32
    algo = 'p';lmax = 4; diffModParamsType = [2,0];
  case 3+32
    algo = 'p';lmax = 6; diffModParamsType = [2,0];
  case 4+32
    algo = 'p';lmax = 8; diffModParamsType = [2,0];
  case 5+32
    algo = 'p';lmax = 10;diffModParamsType = [2,0];
  case 6+32
    algo = 'p';lmax = 12;diffModParamsType = [2,0];
  case 7+32
    algo = 'p';lmax = 14;diffModParamsType = [2,0];
  case 8+32
    algo = 'p';lmax = 16;diffModParamsType = [2,0];
      
    % Small tensor canonical
  case 9+32
    algo = 'p';lmax = 2; diffModParamsType = [1,.11];
  case 10+32
    algo = 'p';lmax = 4; diffModParamsType = [1,.11];
  case 11+32
    algo = 'p';lmax = 6; diffModParamsType = [1,.11];
  case 12+32
    algo = 'p';lmax = 8; diffModParamsType = [1,.11];
  case 13+32
    algo = 'p';lmax = 10;diffModParamsType = [1,.11];
  case 14+32
    algo = 'p';lmax = 12;diffModParamsType = [1,.11];
  case 15+32
    algo = 'p';lmax = 14;diffModParamsType = [1,.11];
  case 16+32
    algo = 'p';lmax = 16;diffModParamsType = [1,.11];

    
  %% Deterministic tractography
  % Beherens (1,0)
  case 1+48
    algo='d';lmax = 2;diffModParamsType = [1,0];
  case 2+48
    algo='d';lmax = 4;diffModParamsType = [1,0];
  case 3+48
    algo='d';lmax = 6;diffModParamsType = [1,0];
  case 4+48
    algo='d';lmax = 8;diffModParamsType = [1,0];
  case 5+48
    algo='d';lmax = 10;diffModParamsType = [1,0];
  case 6+48
    algo='d';lmax = 12;diffModParamsType = [1,0];
  case 7+48
    algo='d';lmax = 14;diffModParamsType = [1,0];
  case 8+48
    algo='d';lmax = 16;diffModParamsType = [1,0];
 
    % Optimal tensor form Corpus Callosum
  case 9+48
    algo='d';lmax = 2; diffModParamsType = [1.8,0.2];
  case 10+48
    algo='d';lmax = 4; diffModParamsType = [1.8,0.2];
  case 11+48
    algo='d';lmax = 6; diffModParamsType = [1.8,0.2];
  case 12+48
    algo='d';lmax = 8; diffModParamsType = [1.8,0.2];
  case 13+48
    algo='d';lmax = 10;diffModParamsType = [1.8,0.2];
  case 14+48
    algo='d';lmax = 12;diffModParamsType = [1.8,0.2];
  case 15+48
    algo='d';lmax = 14;diffModParamsType = [1.8,0.2];
  case 16+48
    algo='d';lmax = 16;diffModParamsType = [1.8,0.2];
   
   % isotropic tensor (1,1)
  case 1+16+48
    algo='d';lmax = 2;diffModParamsType = [1,1];
  case 2+16+48
    algo='d';lmax = 4;diffModParamsType = [1,1];
  case 3+16+48
    algo='d';lmax = 6;diffModParamsType = [1,1];
  case 4+16+48
    algo='d';lmax = 8;diffModParamsType = [1,1];
  case 5+16+48
    algo='d';lmax = 10;diffModParamsType = [1,1];
  case 6+16+48
    algo='d';lmax = 12;diffModParamsType = [1,1];
  case 7+16+48
    algo='d';lmax = 14;diffModParamsType = [1,1];
  case 8+16+48
    algo='d';lmax = 16;diffModParamsType = [1,1];
 
    % Small tensor
  case 9+16+48
    algo='d';lmax = 2; diffModParamsType = [.5,0];
  case 10+16+48
    algo='d';lmax = 4; diffModParamsType = [.5,0];
  case 11+16+48
    algo='d';lmax = 6; diffModParamsType = [.5,0];
  case 12+16+48
    algo='d';lmax = 8; diffModParamsType = [.5,0];
  case 13+16+48
    algo='d';lmax = 10;diffModParamsType = [.5,0];
  case 14+16+48
    algo='d';lmax = 12;diffModParamsType = [.5,0];
  case 15+16+48
    algo='d';lmax = 14;diffModParamsType = [.5,0];
  case 16+16+48
    algo='d';lmax = 16;diffModParamsType = [.5,0];
      
    % Large tensor
  case 1+32+48
    algo='d';lmax = 2; diffModParamsType = [2,0];
  case 2+32+48
    algo='d';lmax = 4; diffModParamsType = [2,0];
  case 3+32+48
    algo='d';lmax = 6; diffModParamsType = [2,0];
  case 4+32+48
    algo='d';lmax = 8; diffModParamsType = [2,0];
  case 5+32+48
    algo='d';lmax = 10;diffModParamsType = [2,0];
  case 6+32+48
    algo='d';lmax = 12;diffModParamsType = [2,0];
  case 7+32+48
    algo='d';lmax = 14;diffModParamsType = [2,0];
  case 8+32+48
    algo='d';lmax = 16;diffModParamsType = [2,0];
      
    % Small tensor canonical with the same ration md/rd (1.8/0.2=1/.11)
  case 9+32+48
    algo='d';lmax = 2; diffModParamsType = [1,.11];
  case 10+32+48
    algo='d';lmax = 4; diffModParamsType = [1,.11];
  case 11+32+48
    algo='d';lmax = 6; diffModParamsType = [1,.11];
  case 12+32+48
    algo='d';lmax = 8; diffModParamsType = [1,.11];
  case 13+32+48
    algo='d';lmax = 10;diffModParamsType = [1,.11];
  case 14+32+48
    algo='d';lmax = 12;diffModParamsType = [1,.11];
  case 15+32+48
    algo='d';lmax = 14;diffModParamsType = [1,.11];
  case 16+32+48
    algo='d';lmax = 16;diffModParamsType = [1,.11];
    
  otherwise
    keyboard
end
  end