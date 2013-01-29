function   fg = mctTrackInRoi(dt,roi,numFibersPerCoord,fgName,algo)
% function fg = mctTrackInRoi(dt,roi,numFibersPerCoord,[fgName],[algo])
%
% This will track every fiber in ROI.
% 
% dt is a structure.
%
% See also: t_fgCompute, mctDiffusionModel, dwiGet, fefgGet
%
% Example:
%  See t_fgCompute
%   [A dSig] = mctDiffusionModel(dwi,fgImg,theseCoords);
%    w       = mctDiffusionModeFit(A,dSig);
%
% Franco
%
% (c) Stanford VISTA Team 2012

% This is the standarddeviation of the gaussian distribution used to
% generate random seeds around each coordinate in coords. Units are those
% in mmPerVoxel.
stdOfSeedsJitter = 0; % In mm

% Distance between steps in the tractography algoithm
opts.stepSizeMm  = 1; % by keeping this set to 1 we sample one fiber per eantry in seeds.
                      % seeds are coords repmat the number of fibers

% Stopping criteria FA<0.2
opts.faThresh    = 0.2;

% Stopping criteria angle between steps >30 degrees
opts.angleThresh = 45;

% Uknown ... (Bob?)
opts.wPuncture   = 0.2;

% Interpolation method. After each step we interpolate the tensor at that
% point. Trilinear interpolation works well.
opts.whichInterp = 1;

% Discard Fibers shorter than 2mm or longer than 250mm
opts.lengthThreshMm = [2 250];

% There are multiple algorithms that can be used.  We prefer STT. See:
% Basser PJ, Pajevic S, Pierpaoli C, Duda J, Aldroubi A. 2000.
% In vivo fiber tractography using DT-MRI data.
% Magnetic Resonance in Medicine 44(4):625-32.
if notDefined('algo'), algo = 1; end
opts.whichAlgorithm = algo;

% Perform tractography
fg = dtiFiberTrack2(dt.dt6, roi.coords', numFibersPerCoord, stdOfSeedsJitter, dt.mmPerVoxel, dt.xformToAcpc, fgName, opts);

% Show the fibers and the ROI
% figure
% for iff = 1:length(fg.fibers)
%   plot3(round(fg.fibers{iff}(1,:)),round(fg.fibers{iff}(2,:)),round(fg.fibers{iff}(3,:)),'.-');
%   hold on
% end
% 
% plot3(roi.coords(:,1),roi.coords(:,2),roi.coords(:,3),'ro','MarkerFaceColor','r','MarkerSize',10)
% plot3(roi.coords(:,1),roi.coords(:,2),roi.coords(:,3),'ks','MarkerFaceColor','k','MarkerSize',10)
% 
% set(gca,'XTick',[-81:2:80],'ZTick',[-81:2:80],'YTick',[-81:2:80])
% grid on


