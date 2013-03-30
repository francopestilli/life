function done = s_ms_check_processes(feFileToLoad,trackingType,lmax,bval,cullType)
%
% s_ms_test_check_processes(feFileToLoad,trackingType,lmax,bval,culleType))
%
% Check the FE structures for a series of connectomes.
% Shows which processes were run on the connectomes, e.g., culling.
%
% This is part of a series of reproducible science scritps to be published
% with the LiFE mansucript.
%
% Franco (c) Stanford Vista Team 2013
if notDefined('trackingType'),trackingType = 't';end
if notDefined('lmax'),        lmax         = 2;end
if notDefined('bval'),        bval         = 2000;end
if notDefined('rep'),         rep          = [1,2,3];end
if notDefined('diffusionModelParams'),   diffusionModelParams=[1,0];end
if notDefined('cullType'),   cullType='culledL2';end
if notDefined('feFileToLoad'), feFileToLoad = 'load from disk';
else rep =1;
end
fprintf('[%s]\n',mfilename)
done = zeros(length(rep),1);

for ilmx = 1:length(lmax)
culledRep = cell(3,1);
bytes     = zeros(3,1);
for irep = 1:length(rep)
    % Get the fe structure
    if ischar(feFileToLoad)
        % Information on the path to the files to load.
        % This is where the inputs will be loaded from
        [feFileToLoad, fname] = ...
            msBuildFeFileName(trackingType,lmax(ilmx),bval,rep(irep), ...
            diffusionModelParams,cullType);
        
        if exist(feFileToLoad,'file')
            %fprintf('Found %s ...\n',feFileToLoad)
            fileInfo = dir(feFileToLoad);
            
            if fileInfo.bytes > 3000
                culledRep{irep} = 'done';
                done(irep) = 1;
            else
                culledRep{irep} = 'in progress';
            end
            bytes(irep) = fileInfo.bytes;
            
        else
            culledRep{irep} = 'not started';
        end
    end
end
fprintf('File:\n%s\nCulling information:\n',fname);
for irep = 1:length(rep)
    fprintf('---> <%i> <%s> <%i> bytes.\n',irep,culledRep{irep},bytes(irep))
end
fprintf('\n')
end

end