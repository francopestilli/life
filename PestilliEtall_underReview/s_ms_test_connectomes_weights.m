function s_ms_test_connectomes_weights(feFileToLoad,trackingType,lmax,bval,rep)
%
% s_ms_test_connectomes_weights(trackingType,lmax,bval,rep)
%
% Test the hypothesis that a ascicle is important for each fascicle int he
% connectome.
%
% The code has not been worked out well enough yet.
%
% This is part of a series of reproducible science scritps to be published
% with the LiFE mansucript.
%
% Franco (c) Stanford Vista Team 2013
if notDefined('trackingType'),trackingType = 'p';end
if notDefined('lmax'),        lmax         = 8;end
if notDefined('bval'),        bval         = 4000;end
if notDefined('rep'),         rep          = [1,2,3];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('cullType'),   cullType='culledL2';end

if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','connectomes_weights');end
% High-resolution Anatomy
t1File = '/azure/scr1/frk/150dirs_b1000_b2000_b4000/150dirs_b2000/t1/t1.nii.gz';

% Check that the connectome were proprocessed before attempting to make a
% plot.
done = s_ms_check_processes([],trackingType,lmax,bval,cullType);

if ~all(done), 
    fprintf('\n[%s] Not all connectomes were proprocessed... not proceeding with plot.\n',mfilename);
return
end

for irep = 1:length(rep)
        % Information on the path to the files to load.
        % This is where the inputs will be loaded from
        [feFileToLoad, fname] = msBuildFeFileName(trackingType,lmax,bval,rep(irep),diffusionModelParams,cullType);

    
    % Get the fe structure
    disp('loading the LiFE structure...')
    if ischar(feFileToLoad)
        fprintf('Loading %s ...\n',feFileToLoad)
        load(feFileToLoad);
    else
        fe  =feFileToLoad;
        clear feFileToLoad;
    end
    weights = feGet(fe,'fiber weights');
    
    % Remove one fiber at the time and compute the variance explained with
    % and withtu that fiber
    for ifib = 1:feGet(fe,'nfibers')
        fasIndex       = false(feGet(fe,'nfibers'),1);
        fasIndex(ifib) = 1;
        
        % Remove all the voxels from the connectome except the ones where the
        % fascicle passes through. Fit the new model.
        [feWithoutFas, feWithFas, con(ifib).c] = feTestFascicle(fe,fasIndex,0);
        
        w.rmse{ifib}     = (feGetRep(feWithFas,   'vox  rmse'));
        wo.rmse{ifib}  = (feGetRep(feWithoutFas,'vox  rmse'));
        
        % Get a probability value for the connection represented by this
        % fascicle
        [p(ifib), EmpiricalDiff(ifib)] = feComputeConnectivity(w.rmse{ifib},wo.rmse{ifib});
        fprintf('[%s]\n  Pc < %2.2f%% \n  Wf:%2.4f \n  Ed:%2.4f\n', ...
            mfilename,p(ifib),weights(ifib),EmpiricalDiff(ifib))

    end
    keyboard
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName)
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end
eval( sprintf('print(%s,  ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'');', num2str(h),figName));

end
