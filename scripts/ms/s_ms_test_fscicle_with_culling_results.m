function s_ms_test_fscicle_with_culling_results(fitType,baseDir)
%
% This function makes some plots of the results of a culling of connectomes
% across iterations.
%
% It is in dev state, it might improve if these plots become useful.
%
% Franco Pestilli (c) Stanford University 2013

close all
clear h

if notDefined('fitType'), fitType = 'L1';end
os  = {sprintf('oTensor%s',fitType),sprintf('oDet%s',fitType),sprintf('oProb%s',fitType)};
fes = {sprintf('feTensor%s',fitType),sprintf('feDet%s',fitType),sprintf('feProb%s',fitType)};

% Load the results from disk:
for ii = 1:length(fes)
    if ~exist(fes{ii},'var')
        disp('loading culling results...')
        load(fullfile('/home/frk/',[fes{ii}(1:end-2),'.mat']))
    end
end

% Data directories:
% used copy: /home/frk
if notDefined('baseDir')
    baseDir = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/results/life_mrtrix_rep1/culling_results_l1_l2/';
end

s_ms_test_connectomes_hypothesis(feProbL2)

keyboard