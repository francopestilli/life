function [strand_collection, isotropic_regions] = mctNfgLoadStrands(strandsDir)
% Load NFG strands into Matlab struct.
%
%   strand_collection = mctNfgLoadStrands(strandsDir)
%
% Example:
%  fg = mctNfgLoadStrands(dirname)
%
%  strandDir is the directory containing the NFG strands (fibers)
%
% by Franco Pestilli
%
% (C) 2012 Stanford VISTA team. 

if notDefined('strandsDir'), strandsDir = '.'; end

strand_counter    = 0;
isotropic_regions = [];

files = dir(strandsDir)'; %'maybe here i need to load onluy the files starting by strand*??
if (size(files) == [1,0])
    error(sprintf('[%s] Could not load any strands from directory %s.', mfilename, strandsDir));
end

% if various files are present in the directory, 
% subselect the strands files only
for file = files   
    delimiters = [strfind(file.name, '_') strfind(file.name, '-') strfind(file.name, '.txt' )];
    if  ~isempty(regexp(file.name,'strand', 'once')) 
        strand_counter = strand_counter + 1;        
        strands(strand_counter) = file;
   
        % this sorting code was take from NFG display_strands.m
        strand_collection{strand_counter, 1} = load([strandsDir filesep file.name]);
        strand_collection{strand_counter, 2} = str2num(file.name(delimiters(1) + 1 :delimiters(2) -1 ));
        strand_collection{strand_counter, 3} = str2num(file.name(delimiters(3) + 2 :delimiters(4) -1 ));
        strand_collection{strand_counter, 4} = str2num(file.name(delimiters(2) + 1 :delimiters(3) -1 ));
    end
    if ~isempty(regexp(file.name,'iso', 'once'))
        isotropic_regions = load([strandsDir filesep file.name]);
    end
end

return