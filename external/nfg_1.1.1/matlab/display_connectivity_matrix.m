function [matrices, not_null_bundle_indices] = display_connectivity_matrix(varargin)

% function [bundle_indices] = display_connectivity_matrix(dirname, varargin)
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


tick_increment = 4;

order_by_column = 1;

rois_info_supplied = 0;
num_matrices = 0;

matrix_names = {};

if (nargin == 0)
	error('No connectivity matrix paths supplied');
else

    
	combined_not_null_rois = [];

	for (varargin_i = 1:nargin)
      
       
        if (strcmp(varargin{varargin_i},'indices'))
        
            order_by_column = 1;
            
        elseif (strcmp(varargin{varargin_i},'cross_section'))
        
            order_by_column = 3;
            
        elseif (strcmp(varargin{varargin_i},'length'))
        
            order_by_column = 4;    
            
        elseif (size(strfind(varargin{varargin_i}, 'rois_info.txt'),1) ~= 0) 
            
            rois_info_supplied = 1;
            
            bundle_info_matrix = load(varargin{varargin_i});
            
		else 

			
            matrix = load(varargin{varargin_i});
			
            num_matrices = num_matrices + 1;
			
            matrix_names{end+1} = varargin{varargin_i};
            
            if ((size(matrix,1)+1) ~= size(matrix,2))
                error(['Error loading matrix from ' path ', number of rows (' num2str(size(matrix,1)) ') does not match number of columns (' num2str(size(matrix,2)) '). ']);
            end

            not_null_rois = [];

            for (row_i = 2:size(matrix,1))

                not_null_elems = find(matrix(row_i,:));

                num_row_elems = size(not_null_elems,2);	

                if (num_row_elems > 0)
                    not_null_rois(end+1) = row_i;
                end
            end

            not_null_rois = sort(not_null_rois);


            matrices{num_matrices} = matrix; %#ok<AGROW>
            combined_not_null_rois = [combined_not_null_rois; not_null_rois']; %#ok<AGROW>
            not_null_bundle_indices{num_matrices} = unique(floor(not_null_rois./2) - 1); %#ok<AGROW>
            

        end
		
    end
		
    

    
    if (rois_info_supplied)
        
        bundle_info_matrix = sortrows(bundle_info_matrix, -order_by_column);
        
        rois_to_be_included = [];
        
        for (bundle_i = 1:size(bundle_info_matrix,1))
            
            rois_to_be_included(end+1) =  bundle_info_matrix(bundle_i,5);
            rois_to_be_included(end+1) =  bundle_info_matrix(bundle_i,6);
            
        end
        
        
        
    else
        
        rois_to_be_included = sort(unique(combined_not_null_rois));
        
    end
	

    figure();


    for (matrix_i = 1:num_matrices)
        
		matrix = matrices{matrix_i};
		
		matrix = matrix(rois_to_be_included,:);
		matrix = matrix(:,[(rois_to_be_included+1), 1, 2]);

        subplot(1,num_matrices,matrix_i);
        
        imagesc(matrix, [0 1]);
        daspect([1,1,1]);
        colorbar;
	
		
        if (rois_info_supplied)
	
			ticks = [1,tick_increment:tick_increment:size(bundle_info_matrix,1)]; 
			
			set(gca, 'YTick', (ticks .* 2));
            set(gca, 'YTickLabel', num2str(bundle_info_matrix(ticks,order_by_column),2));

            set(gca, 'XTick', ([ticks,(ticks(end)+tick_increment)] .* 2));         
            set(gca, 'XTickLabel', num2str([bundle_info_matrix(ticks,order_by_column);''],2));
        end
        
        title(matrix_names(matrix_i));
         
    end


	

	
	
	
	
end

	
	
