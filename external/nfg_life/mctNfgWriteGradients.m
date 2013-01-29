function fid = mctNfgWriteGradients(nfgGradientFile,bvecs,bvals)
%
% function fid = mctnfgWriteGradients(nfgGradientFile,bvecs,bvals)
%
% writes bvecs and bvals into an NFG gradient.txt file
%
% Franco
%
% this needs some more work to figure out how to save the gradients in the order NFG wants them.
% but I am not 100% sure i need this function at all.


% create gradients with the correct orientation
gradients = [bvecs,bvals];
if size(gradients,2) == 4; gradients = gradients';
end

fid = fopen(nfgGradientFile,'wt');
fprintf(fid,'%1.4f %1.4f %1.4f %1.0f\n',gradients)
fclose(fid);
