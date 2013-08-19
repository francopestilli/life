function s_ms_test_connectomes_volume(trackingType,lmax,diffusionModelParams)
%
% Load FE structeres obtained by preprocessing connectomesconstrained within a
% region of interest and within the cortex and makes some basic plot of
% statistics in the connectomes
%
% See also:
%
% Franco (C) 2012 Stanford VISTA team.

if notDefined('trackingType'),trackingType = 'p';end
if notDefined('lmax'),        lmax         = 8;end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = [1,2];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('cullType'),   cullType='culledL2';end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_volume');end
% Check that the connectome were proprocessed before attempting to make a
% plot.
% done = s_ms_check_processes([],trackingType,lmax,bval,cullType);
% if ~all(done), 
%     fprintf('\n[%s] Not all connectomes were proprocessed... not proceeding with plot.\n',mfilename);
%     return
% end

for i_lmax = 1:length(lmax)
    % Loop over Bvalues and tracking reps.
    for i_bval = 1:length(bval)
        for irep = 1:length(rep)
            
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
            
            % Make a plot of volumetric analyse
            [percentVol(irep), figH] = fePlot(fe,'wmvolume');
            for fi = 1:length(figH)
                figName = sprintf('volume_%s',fname);
                
                % Save the figure
                saveFig(figH(fi),fullfile(saveDir,figName))
            end
            close all; drawnow
        end
        
        
        % Make a plot across repetitions:
        % Volume of the ROI
        figName = sprintf('Percent_Volume_%s',fname);
        fh = mrvNewGraphWin(figName);
        h = bar(mean(percentVol),'k');
        hold on
        plot([1,1],[mean(percentVol), mean(percentVol)]+[-std(percentVol), std(percentVol)],'r-','LineWidth',10)
        set(gca,'ylim',[60 110],'xlim',[.5 1.5])
        ylabel('Percent white-matter volume')
        title(sprintf('Percent volume: %2.3f (+/-%2.3f)',mean(percentVol),std(percentVol)))
        set(gca,  'ytick',[60 80 100 120],  'box','off','tickDir','out')
        % Save the figure
        saveFig(fh,fullfile(saveDir,figName))
    end
end

end % Main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)

fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end
eval( sprintf('print(%s,  ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'');', num2str(h),figName));

end