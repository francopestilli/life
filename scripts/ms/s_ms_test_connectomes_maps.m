function s_ms_test_connectomes_maps(trackingType,dataType,lmax,diffusionModelParams,recompute)
%
% Load FE structeres obtained by preprocessing connectomesconstrained within a
% region of interest and within the cortex and makes some basic plot of
% statistics in the connectomes
%
% See also:
%
% Franco (C) 2012 Stanford VISTA team.

if notDefined('trackingType'),trackingType = 't';end
if notDefined('lmax'),        lmax         = 6;end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = [1,2,3];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('cullType'),   cullType='culledL2';end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_hists_maps');end
% Check that the connectome were proprocessed before attempting to make a
% plot.
done = s_ms_check_processes([],trackingType,lmax,bval,cullType);
if ~all(done),
    fprintf('\n[%s] Not all connectomes were proprocessed... not proceeding with plot.\n',mfilename);
    return
end
for irep = 1:length(rep)
    for i_lmax = 1:length(lmax)
        % Loop over Bvalues and tracking reps.
        for i_bval = 1:length(bval)
            
            % File to load
            % This is where the inputs will be loaded from
            [feFileToLoad, fname] = ...
                msBuildFeFileName(trackingType,lmax(i_lmax),bval(i_bval),rep(irep), ...
                diffusionModelParams,cullType);
            
            if (exist(feFileToLoad,'file') == 2)
                fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)
                load(feFileToLoad);
            else
                fprintf('[%s] FE file NOT found: \n%s\n ======================================== \n\n\n',mfilename,feFileToLoad)
                keyboard
            end
            
            if irep == 1
            % Make a plot of the maps
            slice = 26:35;
            for is = 1:length(slice)
                [~, figH] = fePlot(fe,'fiberdensitymap',slice(is));
                for fi = 1:length(figH)
                    figN = get(figH(fi),'name');
                    figName = sprintf('%s_slc%i_lmx%s',figN(~isspace(figN)),slice(is),num2str(lmaxs(i_lmax)));
                    
                    % Save the figure
                    saveFig(figH(fi),fullfile(saveDir,figName))
                end
            end
            end
            
            % Make a plot of the hists
            [~, figH] = fePlot(fe,'fiberdensityhist');
            for fi = 1:length(figH)
                figN = get(figH(fi),'name');
                figName = sprintf('%s_lmx%s',figN(~isspace(figN)),num2str(lmaxs(i_lmax)));
                
                % Save the figure
                saveFig(figH(fi),fullfile(saveDir,figName))
            end
            close all; drawnow
        end % irep
    end   % i_bval
end

end % Main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end

printCommand = sprintf('print(%s, ''-djpeg90'', ''-opengl'', ''%s'')', num2str(h),figName);
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

% do the printing here:
eval(printCommand);

% save a wiki-compatible file
eval(sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName));

end