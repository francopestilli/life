function [fe, cull] = feConnectomeCull(fe,display)
%
% fe = feConnectomeCull(fe)
%
% Fits a connectome over and over by removing the fibers that make no
% contribution to the diffusion signal at each iteration.
%
% Copyright Franco Pestilli, Stanford University 2014

if notDefined('display'); display.cull = false;end

% Initialize a structure containing the information abotu the culling
cull.name = feGet(fe,'name');
cull.iter = 1;
cull.minwval = eps; % (here we could have used 0 as thresholds, we use eps instead 
                    % in hopes to seepd things up)

% Extract the connectome model and the signa for fitting.
M    = feGet(fe,'mfiber');
dSig = feGet(fe,'dsigdemeaned');
fe   = feSet(fe,'fit',feFitModel(M,dSig,'bbnnls'));

% Count the number of fbers that make no contribution to the diffusion
% prediction (here we could have used 0 as thresholds, we use eps instead 
% in hopes to seepd things up)
sw   = sum(feGet(fe,'fiber weights') <= cull.minwval);

% Update the culling structure given the results of the firts fit.
cull.num2delete(cull.iter) = sw;
cull.numtotal(cull.iter)   = size(M,2);
cull.rmse(cull.iter)       = mean(feGet(fe,'vox rmse'));
cull.rmse(cull.iter)       = mean(feGetRep(fe,'vox rmse'));
cull.rrmse(cull.iter)      = mean(feGetRep(fe,'vox rmse ratio'));

% Cull the connectome. We stop when no fibers have a contribution to the
% diffusion signal less than eps.
while (sw ~= 0)
    fprintf('Number of total %i and zero-weight %i fibers | culling...\n', ...
        cull.numtotal(cull.iter),cull.num2delete(cull.iter))
    M    = feGet(fe,'mfiber');
    dSig = feGet(fe,'dsigdemeaned');
    fe   = feSet(fe,'fit',feFitModel(M,dSig,'bbnnls'));
    sw   = sum(feGet(fe,'fiber weights') <= cull.minwval);
    fe   = feConnectomeReduceFibers(fe, find((feGet(fe,'fiber weights') > cull.minwval)));
    cull.iter = cull.iter + 1;
    cull.num2delete(cull.iter) = sw;
    cull.numtotal(cull.iter)   = size(M,2);
    cull.rmse(cull.iter)       = mean(feGet(fe,'vox rmse'));
    cull.rmse(cull.iter)       = mean(feGetRep(fe,'vox rmse'));
    cull.rrmse(cull.iter)      = mean(feGetRep(fe,'vox rmse ratio'));
end

% If display cull process
if display.cull
keyboard    
end

end