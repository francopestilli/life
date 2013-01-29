function display_ext_stats(collection_dir, varargin)

%	function display_ext_stats('dirname', 'stat_flag', segment_range)
%
%	args: 
%            dirname - name of the input directory (the output directory of a 'mri_sim' call with 'save_ext_stats' set to 1)
%			 varargin - Limit the plots displayed to either
%	displays: A slice by slice representation of the plotted subvoxel orientations.  
%
%
%	NB: To generate the extended statistics, set the input parameter 'save_ext_stats' to 1 and run 'mri_sim' on a strand collection. 
%         
%
%	display_ext_stats
%	Numerical fibre Generator
%
%
%	Created by Tom Close on 17/07/07.
%   Copyright 2008 Tom Close.
%   Distributed under the GNU General Public Licence.
%
% 
% 
%   This file is part of 'Numerical Fibre Generator'.
% 
%   'Numerical Fibre Generator' is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   'Numerical Fibre Generator' is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with 'Numerical Fibre Generator'.  If not, see <http://www.gnu.org/licenses/>
% 
%

 	if (nargin < 1) 
 	   error('Incorrect number of arguments:\n Usage:\n    arg[0] = dirname\n    arg[2]= stat_flag\n    arg[3]= segment_range');
	elseif (nargin == 1)
		stats_type_list = {'segment_length', 'segment_angle', 'segment_overlap_fraction', 'segment_radius_curv'};
		num_types = 4;
		
	else
		stats_type_list = varargin;
		num_types = nargin - 1;
		
	end


	
	
	
%	new_strand_markers = load_ext_stats(collection_dir, 'new_strand_markers.int');
	
%	new_strand_indices = find(new_strand_markers);
%	strand_labels = [1:size(new_strand_indices)];
	
	for (type_i = 1:num_types) 
		
		if (~ischar(stats_type_list{type_i}))
			error(['Arg ' num2str(type_i+1) ' is not a valid string (it needs to be the name of a statisical output file i.e. "segment_length", "segment_angle", "segment_overlap_fraction" or "segment_radius_curv" )']);
		end	
		
		stats_type = stats_type_list{type_i};
		
		stats = load_ext_stats(collection_dir, [stats_type, '.double']); 
		
		stats = sort(stats);
		
		second_percentile_index = ceil(size(stats,1) * 0.02);
		ninety_eighth_percentile_index = floor(size(stats,1) * 0.98);
		
		second_percentile = stats(second_percentile_index);
		ninety_eighth_percentile = stats(ninety_eighth_percentile_index);	
		
		increment = (ninety_eighth_percentile - second_percentile)/200;
		
		bins = [max( [0 min(stats)]):increment:min( [(ninety_eighth_percentile + 800 * increment) max(stats)])];
		      
 		h = figure();
 		set(h,'Units','pixels') 
 		set(h, 'Position', [750 300 1200 740]);
     
    	hist(stats, bins)
       
         if (strcmp(stats_type, 'segment_radius_curv'))          
 
             xlabel('Approx. radius of curvature*  ', 'FontSize', 11, 'FontWeight', 'bold');
             set(gca, 'Xlim', [0 2.0]); 
             
 		elseif (strcmp(stats_type, 'segment_length'))
    
            xlabel('Control point interval length* ', 'FontSize', 11, 'FontWeight', 'bold');
            
        elseif (strcmp(stats_type, 'segment_angle'))
    
            xlabel('Angular deviation at control points (degrees) ', 'FontSize', 11, 'FontWeight', 'bold');
            
        elseif (strcmp(stats_type, 'segment_overlap_fraction'))
    
            xlabel('Fraction of segment that is overlapping ' , 'FontSize', 11, 'FontWeight', 'bold');
        end        
              
		
	end
	
end


function stats = load_ext_stats(dirname, stats_filename)

	stats_dir = [dirname '/stats/'];

	if (~isdir(stats_dir)) 
		error(['No statistics directory found for ' dirname]);
		
	end


	
	[stats_type, file_ext] = strtok(stats_filename, '.');
	
	data_type = file_ext(1,(2:size(file_ext,2)));
	
	filename = [stats_dir stats_filename];
	
	file = fopen(filename);

	if (file == -1)
		error([ 'Could not file ' filename '.  Skipping...\n' ]);
	end
		
	stats = fread(file, data_type);

	fclose(file);



end

