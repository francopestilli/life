function [dwi dwiNifti dwiNiftiFile] = mctNfgBuildDwi(dwiDir)
% Load a series of analyze (3D) files and builds and saves a (4D) nifti.
%
% function [dwiNifti dwi] = mctNfgBuildDwi(dwiDir)
%
% Example:
%
%   fg = mctNfgBuildDwi(dwiDir);
%
% See also: mctNfgSimulate.m, mctNfgLoadStrands.m, s_mct_simulated.m
%
% by Franco Pestilli
%
% (C) 2012 Stanford VISTA team. 

% Load the DWI data.
if ~exist(fullfile(dwiDir,'dwi_all.nii.gz'),'file')
    % get all the file names
    files = dir(dwiDir)';
    
    % did we find any?
    if (size(files) == [1,0])
        error(sprintf('[%s] Could not find any DWI image in directory %s.', mfilename, strandsDir));
    end
        
    % (2) load one file, to use as a template for the nifti. clear data change
    % the filename
    done = 0;
    for file = files
        if  ~isempty(regexp(fullfile(dwiDir,file.name),'.img', 'once')) && ~done
            dwiNifti = niftiRead(fullfile(dwiDir,file.name));
            dwiNifti.data  = [];
            dwiNifti.fname = sprintf('dwi_all.nii.gz');
            dwiNifti.ndim  = 4;
            done = 1;
        end
    end
    
    
    % either load the individual ANLYZE files created by NFG (and contatenate them)
    for ii = 1:length(files)
        if  ~isempty(regexp(fullfile(dwiDir,files(ii).name),'.img', 'once'))
            temp = niftiRead(fullfile(dwiDir,files(ii).name));
            dwiNifti.data = cat(4, dwiNifti.data, temp.data);
            % fprintf('[%s] loading: %s tempVal = %2.5f dataVal = %2.5f\n',mfilename,files(ii).name,temp.data(1,1,1),squeeze(dwiNifti.data(1,1,1,end)))
        end
    end
    dwiNifti.dim = [dwiNifti.dim size(dwiNifti.data,4)];
    
    % (4) save the nifti file in a good place
    dwiNiftiFile = fullfile(dwiDir,dwiNifti.fname);
    niftiWrite(dwiNifti,dwiNiftiFile);
    
else
    % or load the load the nifti file created by mct if this is not the first time we use mctNfgBuildDwi on this simulation
    dwiNifti = niftiRead(fullfile(dwiDir,'dwi_all.nii.gz'));
end

% (5) load the bvecs and bvals and save them in a dwi structurre.
% the location of the bvecs is assumed to be ../parameters
% this folder is saved by default when using the bash shell script 'run'
bvecsFile  = 'grad_directions.txt'; 
delimeters = [strfind(bvecsFile, '_') strfind(bvecsFile, '-') strfind(bvecsFile, '.txt' )];

bvecs = load(fullfile(dwiDir,'parameters',bvecsFile));
bvals = bvecs(:,4)/1000;
bvecs = bvecs(:,1:3);

% (5)  build a dwi file using the nifit file and the bvals, bvecs
dwi    = dwiCreate('name','dwi-struct','nifti',dwiNifti,'bvecs',bvecs,'bvals',bvals);


