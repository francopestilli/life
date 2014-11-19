function fe = feSet(fe,param,val,varargin)
% Set fascicle evaluation parameters.
%
%   fe = feSet(fe,param,val,varargin)
%
%----------
% feSet(fe,'bvecsindices');
%----------   
% Set the a second diffusion signla measuremnt (repeat).
% Image_vals is the image array with the dwi signal taken from all the
% voxels in the new dataset the fe.roi
% fe = feSet(fe,'diffusion signal repeat',image_vals);
%----------
% Set the the s0 (no diffusion direction measurement) for a second
% diffusion signla measuremnt (repeat).
% Image_vals is the image array with the dwi signal taken from all the
% voxels in the new dataset the fe.roi
% fe = feSet(fe,'b0signalrepeat',image_vals);
%----------
%
%
% Copyright (2013-2015), Franco Pestilli, Stanford University, pestillifranco@gmail.com.

% Check for input parameters
if notDefined('fe'),    error('fe structure required'); end
if notDefined('param'), error('param required'); end
if ~exist('val','var'), error('Value required'); end

% Squeeze out spaces and force lower case
param = mrvParamFormat(param);

%%
switch param
  % Book-keeping
  case 'name'
    fe.name  = val;        % This one's name
  case 'type'
    fe.type  = 'faseval'; % Always fascicle evaluation
  case 'savedir'
    fe.path.savedir  = val; % Always fascicle evaluation
  
    % Set top level structure, not just single slot
  case 'life'
    fe  = val;  % Structure of parameters and results from LIFE analysis
    
  case 'fgfromacpc'
    % Fiber group representing the candidate connectome.
    % Everything is in IMAGE coordinates (not in ACPC) in the FE structure. 
    % Here we make this assumption.

    % If it is a string we assume thisis a file name for a fiber group.
    % We load the fiber group from file.
    if ischar(val) && exist(val,'file') 
       [p,n]  = fileparts(val);
       val    = fgRead(val); 

    elseif isstruct(val) && isfield(val,'fibers')
    % We do nothing beside assuming that the fiber group 
    % is passed in ACPC coordinates as default in mrDiffusion
   
    else error('[%s] incorrect val parameter %s.',mfilename, val);
    end
    
    % If a fiber group structure was passed in instead
    % we assume that it was passed in in AC-PC coordinates 
    % (the default in mrDiffusion) we trasform  in IMAGE coordinates
    val = dtiXformFiberCoords(val, fe.xform.acpc2img,'img');
      
    % Finally we set the fiber group in the FE structure in IMAGE coordinates
    fe  = feSet(fe,'fg',val);
    
    % Must clear the fields associated with the fiber group every time a new FG is set.
    fe = feSet(fe,'voxel 2 fiber node pairs',[]);
    
  case {'fgimg', 'fg'}
    % Fiber group candidate fascicles, Connectome.
    % Everything is in img coordinates in LiFE
    %
    fe.fg  = val;
    % Must clear when we change the fg
    fe = feSet(fe,'voxel 2 fiber node pairs',[]);

  case {'fgtensors','tensors','fgq','q'}
    % fe = feSet(fe, 'tensors', tensors);    - set the passed tensors
    fe.fibers.tensors = val;

  case 'roi'
    % Cell array of regions of interest where we evaluate
    % Must clear v2fnp
    fe.roi   = val;
    fe       = feSet(fe,'v2fnp',[]);
    
  case {'roifromfg','fgroi','roifg'}
    name   = sprintf('roi_%s', fe.fg.name);
    randColor = rand(1,3);
    fe.roi = dtiNewRoi(name,randColor,fefgGet(feGet(fe,'fg img'),'unique image coords'));
    
  case 'xform'
    fe.xform = val;  % Transforms between coords for fg, roi, and dwi data
     
  %% Diffusion data related parameters
  case {'bvecs','diffusionbvecs'}
    % feSet(fe,'bvecs');
    fe.bvecs = val;
  case {'bvecsindices','diffusionimagesindicesindwivolume'}
    % feSet(fe,'bvecsindices');
    fe.bvecsindices = val;
  case {'bvals','diffusionbvals'}
    fe.bvals = val;    
  case {'diffusionsignalimage','dsi', 'diffusion_signal_img'}
    fe.diffusion_signal_img = val;
  case {'b0signalimage','b0img', 'diffusion_s0_im','s0image'}
    fe.diffusion_S0_img = val;
  case {'usedvoxels'}
    fe.usedVoxels = val;
  case {'modeltensor'}
    fe.modelTensor = val;
  case {'roivoxels','roicoords'}
    % What space?  What form for the coords?
    % Always in IMG coords in LiFE.
    fe.roi.coords = val;
  
    %% The LiFE model
  case 'mfiber'
    fe.Mfiber = val;             % Fiber portion of M matrix
  case {'measuredsignalfull', 'dsigmeasured'}      % Measured signal in ROI
    fe.dSig  = val;
  case 'fit'
    fe.fit = val;
  case 'voxfit'
    fe.voxfit = val;
  case 'xvalfit'
    fe.xvalfit = val;
  
    %% Connectome fibers information.
  case {'numberofuniquefibersineachvoxel','uniquefibersnum','numberofuniquefibers','numuniquef'}
    fe.fibers.unique.num = val;
  case {'indextouniquefibersineachvoxel','uniquefibersindex','uniqueindex','indexesofuniquefibers','indexuniquef','uniquefibers'}
    fe.fibers.unique.index = val;
  case {'numberoftotalfibersineachvoxel','totalfibernmber','fibersnum','numberoffibers','numf','numfibers'}
    fe.fibers.total.num = val;
  case {'indexoftotalfibersineachvoxel','totalfiberindex','fibersbyvox','fibersinvox'}
    fe.fibers.total.index = val;
  case {'voxel2fibernodepairs','v2fnp'}
    % This has to be cleared whenever we change fg or roi
    fe.voxel2FNpair = val;
    % Spatial coordinate transforms for voxels and fg to coordinate frames
  case {'xformimg2acpc','img2acpc','img2acpcxform'}
    fe.xform.img2acpc = val;
  case {'xformacpc2img','acpc2img','acpc2imgxform'}
    fe.xform.acpc2img = val;
  case {'size','imgsize','volumesize','dims','dim'}
    fe.imagedim = val;
    
    %% Diffusion data reapeted measure parameters
  case 'dwirepeatfile'
    fe.path.dwifilerep = val;  % Diffusion weighted file used for testing results
  case {'diffusionsignalimagerepeat'}
    % Set the a second diffusion signla measuremnt (repeat).
    %
    % Image_vals is the image array with the dwi signal taken from all the
    % voxels in the new dataset the fe.roi
    %
    % fe = feSet(fe,'diffusion signal repeat',image_vals);
    fe.rep.diffusion_signal_img = val;
  case {'s0imagerepeat'}
    % Set the the s0 (no diffusion direction measurement) for a second
    % diffusion signla measuremnt (repeat).
    %
    % Image_vals is the image array with the dwi signal taken from all the
    % voxels in the new dataset the fe.roi
    %
    % fe = feSet(fe,'b0signalrepeat',image_vals);
    fe.rep.diffusion_S0_img = val;
  case {'bvecsrepeat','diffusionbvecsrepeat'}
    % feSet(fe,'bvecsrepeat');
    fe.rep.bvecs = val;
  case {'bvecsindicesrepeat','diffusionimagesindicesindwivolumerepeat'}
    % feSet(fe,'bvecsindicesrepeat');
    fe.rep.bvecsindices = val;
  case {'bvalsrepeat','diffusionbvalsrepeat'}
    % fe = feSet(fe,'bvalsrepeat')
    fe.rep.bvals = val;
  case {'imgsizerepeat'}
    fe.rep.imagedim = val;
    
  case {'anatomyfile'}
    fe.path.anatomy = val;
  case 'dwifile'
    fe.path.dwifile = val;  % Diffusion weighted file used for testing results
  case 'dtfile'
    fe.path.dtfile = val;  % Diffusion weighted file used for testing results

  otherwise
    error('Unknown parameter %s\n',param);
end

end
