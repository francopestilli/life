function record_movie(algo,fibers)
% function record_movie(algo,fibers)
%
% Make movies of various fiber gorups in a brain.
%
% Franco

if notDefined('algo'), algo = 'STT';end % STT or TEND
if notDefined('fibers'),fibers= 'm';end

% Build a name for the fibers I want to load
switch fibers
  case {'mori','mori groups','m'}
    switch algo
      case {'STT','TEND'}
        fibers = sprintf('MoriGroups_clean_D5_L4_%s.mat',algo);
      case {'MRTRIX'}
        fibers = sprintf('MoriGroups_clean_D4_L3_%s.mat',algo);
      otherwise
        keyboard
    end

    % fasclicles to plot
    fasc = {[1 3 5 7 9 11 13 15 17 19],[2 4 6 8 10 12 14 16 18 20]}; % left (odds) and right (even) fascicles
    %fasc = {[1],[8]}; % left (odds) and right (even) fascicles
    sides = [1 2];

  case {'whole','brain','b'}
    fibers = sprintf('WholeBrainFG_%s.mat',algo);
    fasc = {1};
    sides = 1;
    
  otherwise
    keyboard
end

fR = [1:90];
fL = [[65:-3:-99], fliplr([65:-3:-99])];
frames  = {fL,fR};
V = {[-15, 25],[15, 25]};

% cd to the data dir
cd /biac2/wandell6/data/arokem/ModelFits/FP20120420/;
saveDir = '/biac2/wandell6/data/arokem/ModelFits/FP20120420/LiFE/mvs/';
% load some fiber groups
disp('Loading fiber groups...')
load(sprintf('150dirs_b1000_1/fibers/whole_brain_%s/%s',algo,fibers))

% load a t1
disp('Loading T1')
nifti = readFileNifti('t1/t1.nii.gz');

for side = sides % left and right fascicles
  
  % plot a fascicle
  for fi = 1:length(fasc{side})
    n = fasc{side}(fi);
    
    % Prepare the new video:
    vidObj = VideoWriter(fullfile(saveDir,sprintf('%s_%s',fg(n).name,algo)));
    vidObj.FrameRate = 15;
    vidObj.Quality = 100;
    open(vidObj);
    
    fh = mctNfgDisplayStrands(fg(n),[.1 .7 .81]);
    set(fh,'Position',[0 0 .25 .75]);
    hold on
    hm = mctDisplayBrainSlice(nifti,[0 0 -7]);
    axis off
    axis vis3d
    set(gcf,'Color',[0 0 0]);
    set(gca,'LooseInset',get(gca,'TightInset'))
    view(0,90)
    fprintf('\n\n\n[%s] Making movie: %s\n',mfilename, fullfile(saveDir,sprintf('%s_%s',fg(n).name,algo)));
    % add a brain slice and capture the frame
    for i =1:120
      %       h = mctDisplayBrainSlice(nifti,[0 frames{side}(i) 0],[],[],1);
      camorbit(1,0,'data',[.75 1 .75])
      drawnow
      M(i) = getframe(gcf);
      % delete(h)
      drawnow
    end
    fprintf('\n\n\n[%s] Saving movie: %s\n\n\n',mfilename, fullfile(saveDir,sprintf('%s_%s',fg(n).name,algo)));
    
    % erite video and clear memory
    writeVideo(vidObj,M)
    close(vidObj);
    close all; drawnow
    clear M w p fh fh2 hm
    
  end
end
