function [fascicles classification] = feAfqRemoveFascicleOutliers(fascicles,classification,opts)
%
% Removes the outliers in a group of fascicles.
%
% [fascicles classification] = feAfqRemoveFascicleOutliers(fascicles,classification,opts)
%
% INPUTS:
%    fascicles      - a set of fascicles created using feAfqSegment.m
%    classification - a classification strucutre created using
%                     feAfqSegment.m
%    opts           - default options for cleaning the outliers in the
%                     fascicle.
%
% OUTPUTS:
%    fascicles      - a set of fascicles as created by feAfqSegment.m
%                     but withut Outliers.
%    classification - a classification strucutre as created using
%                     feAfqSegment.m, but without outliers
%
% NOTEs: it requres AFQ to be on path.
%
% Franco (c) Stanford Vista Tema 2012

% Parameters for the outliers removal process.
if notDefined('opts')
  opts.stdCutOff   = 3.5;   % Standard deviation fo the 3D gaussian distribution used to
  % represent the fascicle when removing outliers. 3.5 z-scores
  opts.maxLen      = 20;  % Max lenght of fibers to be accepted in cm
  opts.maxNumNodes = 100; % This is used only during the computations does not actually change the nodes
end

% Remove fibers outliers, this creates tighter fascicles.
keep = cell(length(fascicles),1);
fprintf('[%s] Removing outliers from fascicles...\n',mfilename)
for in = 1:length(fascicles)
  fprintf('[%s] %i/%i %s.\n',mfilename,in,length(fascicles),fascicles(in).name)
  [fascicles(in), keep{in}] = AFQ_removeFiberOutliers(fascicles(in),opts.stdCutOff,opts.maxLen,opts.maxNumNodes);
  
  % Find the indices to the current fascicle inside the vector of indices
  % for the whol-brain connectome (fg)
   thisFasIndices = find((classification.index == in));

  % Now remove the fibers that we do not want to keep out of the fascicles.
  classification.index(thisFasIndices(~keep{in})) = 0;
end

return

% The following is a test to show that we are addressing the original
% whole-brain conenctome correctly.
%
% Uncomment the following lines to make sure that we can address the
% facicles with the indices inside classification all the way back to the
% orignal whole-brain connectome.
%
% fas    = fascicles;
% for in = 1:length(fas)
%   fas(in).fibers = {};
%   fas(in).fibers = fg.fibers(classification.index==in);
% end
% Plot the same fascicle fpbtained from the cleaned fibers and addresed
% directly in the whole-brain connectome (fg)
% figure(100);clf
% subplot(2,1,1);plot3(fas(in).fibers{1}(1,:),fas(in).fibers{1}(2,:),fas(in).fibers{1}(3,:));
% hold on
% for ii = 1:length(fas(in).fibers)
%   subplot(2,1,1);plot3(fas(in).fibers{ii}(1,:),fas(in).fibers{ii}(2,:),fas(in).fibers{ii}(3,:));
% end
% title(sprintf('Fascicle addressed directly via indices\ninside the whole-brain connectome'))
% subplot(2,1,2);plot3(fascicles(in).fibers{1}(1,:),fascicles(in).fibers{1}(2,:),fascicles(in).fibers{1}(3,:));
% hold on
% for ii = 1:length(fascicles(in).fibers)
%   subplot(2,1,2);plot3(fascicles(in).fibers{ii}(1,:),fascicles(in).fibers{ii}(2,:),fascicles(in).fibers{ii}(3,:))
% end
% title('Fascicle generated from AFQ remove outliers')
