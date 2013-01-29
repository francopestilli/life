function [combined_mask, bundle_colours, not_null_bundle_indices] = display_rois(dirname, varargin)

%  function [combined_mask, bundle_colours] = display_rois(dirname, varargin)
% 
%   display_rois.m
%   Numerical Fibre Generator
% 
%   Created by Tom Close on 19/02/08.
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
% 

	global bundle_colours;

    if (nargin >= 2)
		bundle_colours = varargin{1};
	else
		bundle_colours = [];
    end


    if (nargin >=3)
       dim = varargin{2};
       
    else
       dim = 3;
    end
    
    if (nargin == 4)
        slice_num = varargin{3};
	else
		slice_num = 0;
	end
    

    [combined_mask, bundle_info] = load_combined_mask(dirname);
	
	roi_ids = unique(combined_mask);
	
	bundle_indices = floor(roi_ids(find(roi_ids > 1))/2) - 1; %#ok<FNDSB>
	
	bundle_indices = sort(unique(bundle_indices));
	
	not_null_bundle_indices = [];
	
	
	for (bundle_i = 0:1:bundle_indices(end))
		
		start_roi_id = 2 * (bundle_i + 1);
		end_roi_id = 2 * (bundle_i + 1) + 1;
		
		start_mask_size = size(find( combined_mask == start_roi_id),1);
		end_mask_size = size(find( combined_mask == end_roi_id),1);

		fprintf(['Bundle ' num2str(bundle_i) '\tsize: start = ' num2str(start_mask_size) '\tend = ' num2str(end_mask_size) '\n']);
		
		if (start_mask_size > 0 || end_mask_size > 0) 
			not_null_bundle_indices = [not_null_bundle_indices; bundle_i]; %#ok<AGROW>
		end
	end
	
% 	if (bundle_indices(1) == -1)
% 		bundle_indices = bundle_indices(2:size(bundle_indices,1));
% 	end
% 	
% 	for bundle_i = 1:size(masks,1)
% 
% 		bundle_indices = [bundle_indices; masks{bundle_i, 2}]; %#ok<AGROW>
% 
% 	end

	not_null_bundle_indices
	
	
	bundle_colours_size = size(bundle_colours,1);
	max_bundle_index = max(bundle_indices);


	if (bundle_colours_size == 0)
		bundle_colours = rand([(max_bundle_index+1), 3]);
	elseif (bundle_colours_size < (max_bundle_index+1))
		error(['Bundle colours matrix does not contain enough colours (' num2str(bundle_colours_size) ') for the loaded rois (' num2str(max_bundle_index) ')']);
	end

	colour_mask = zeros([size(combined_mask) 3]);

	if (nargin < 2)
		display_colour_key(bundle_colours, bundle_indices);
	end
	
	for (x = 1:size(combined_mask,1))
		for (y = 1:size(combined_mask,1))
			for (z = 1:size(combined_mask,1))
				
				roi_id = combined_mask(x,y,z);
				
						
				if (roi_id == 0)			% If the given voxel is empty
					
					colour_code = [0; 0; 0];
					
				elseif (roi_id == 1)		% If the given voxel contains the retainer mask
				
					colour_code = [1;1;1];
					
				else					% Otherwise

					colour_code = squeeze(bundle_colours(floor(roi_id/2),:));
					
				end
				
				colour_mask(x,y,z,:) = colour_code;

			end
		end
	end

	fig = figure();
	
	if (slice_num > 0)
		image(squeeze(colour_mask(slice_num,:,:,:)));  
	else
		for (i = [1:size(colour_mask,dim)])
			figure(fig)
			if (dim == 1) 
				image(squeeze(colour_mask(i,:,:,:)));                  
			elseif (dim == 2)
				image(squeeze(colour_mask(:,i,:,:)));
			else 
				image(squeeze(colour_mask(:,:,i,:)));                    
			end
			title(['Slice ' num2str(i)]);
			pause;
		end
		close(fig);
	end
       

	
end


function [combined_mask, bundle_info] = load_combined_mask(dirname)
% function combined_mask = load_combined_mask(dirname)
%
% loads the combined mask from the directory

	files = dir(dirname)'; %'
	
	filename = [dirname '/rois.img'];

	if (size(files) == [1,0]) 
		error(['Could not load any files from directory ' dirname ]);
	end

    if (size(dir(filename),1) == 0)
		error(['Could not find ' filename]);
	end
				               
	f = fopen(filename, 'r');

	if (f == -1)
	   error(['Could not open file ' filename]); 
	end                


	mask = fread(f, 'int32=>double');

	fclose(f);

	mask_size = size(mask,1);

	dim_size = nthroot(mask_size, 3);

	mask(:) = mask(mask_size:-1:1);

	combined_mask = reshape(mask, dim_size, dim_size, dim_size);
	
	info_filename = [dirname '/rois_info.txt'];	
	
	if (size(dir(info_filename),1) == 0)
		error(['Could not find ' info_filename]);
	end
		
	bundle_info = load(info_filename);
	

end


function masks = load_separate_masks(dirname)
% function masks = load_separate_masks(dirname)
%
% loads the separate masks from the directory (not used currently)

	num_masks = 0;

	files = dir(dirname)'; %'

	if (size(files) == [1,0]) 
		error(['Could not load any files from directory ' dirname ]);
	end
	
	for file = files

		if (~file.isdir)

			delimeters = [strfind(file.name, '_') strfind(file.name, '-') strfind(file.name, '.img' )];
          
            
			if (length(delimeters) == 3 && strmatch('mask', file.name))
				num_masks = num_masks + 1;

				f = fopen([dirname '/' file.name], 'r');
				
                if (f == -1)
                   error(['Could not open file ' dirname '/' file.name]); 
                end
                
				mask = fread(f, 'uint8=>double');
							
                fclose(f);
                
				dim_size = nthroot(size(mask,1), 3);
				
    			masks{num_masks, 1} = reshape(mask, dim_size, dim_size, dim_size);
				masks{num_masks, 2} = str2num(file.name(delimeters(1) + 1 :delimeters(2) -1 ));
				masks{num_masks, 3} = str2num(file.name(delimeters(2) + 1 :delimeters(3) -1 ));

			end

		end

	end
	
	disp(['Loaded ' num2str(num_masks) ' from directory ' dirname ]);
end