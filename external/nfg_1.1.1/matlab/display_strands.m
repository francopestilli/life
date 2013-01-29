function [strand_collection, bundle_colours] = display_strands(varargin)

%	function [strand_collection, bundle_colours] = display_strands(
%                                                           dirname or strand_struct,
%                                                           [bundle_colours],
%                                                           [bundle or strands to_plot],
%                                                           ['bundle' or
%                                                           'strand' (which to plot)],
%                                                            ['3D' or 'line' (for 3D rendered or line plots)],
%                                                           [tube_divs]);
%
%   returns:
%         strand_collection - the loaded strand collection (number of strands X 4 cell array)
%		  bundle_colours - the bundle colours used to display each bundle
%		  for use in subsequent plots
%		  (which are randomly generated if left empty)
%	args: 
%
%   Usage:	
%   
%     arg[1]= Either the directory name of the stored collection or a cell structure containing the
%             loaded strands from a previous call (return variable 1).
% 		arg[2]= (optional) A Nx3 matrix containing the colours in which the different bundles will be
% 		        displayed.
% 		arg[4]= (optional) The indices of the strands/bundles to be displayed.
% 		arg[5]= (optional) Whether the supplied indices correspond to 'strand' or 'bundle' (default 
% 		        'bundle')
% 		arg[6]= (optional) The number of divisions around the diameter of the strand. The default is
%             6, which corresponds to a hexagonal cross-section.
% 
%   Output:	
%   
%     return[1]= A 4xN cell structure containing the strand data (to be used as arg[1] in subsequent
%                calls of the function)
% 		return[2]= A Nx3 matrix containing the colours of the strands displayed (to be used as arg[3]
% 		           in subsequent calls of the function to keep colours consistent)
%
%         
%
%	display_strands
%	Numerical fibre Generator
%
%
%	Created by Tom Close on 17/07/07 (based on the script 'tubeplot' by Anders Sandberg, asa@nada.kth.se, 2005).
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
	
	global bundle_colours;  %bundle colours is made global in case you forget to save the current colours and want to retrieve them.


	if (nargin < 1 || nargin > 6) 
	   error(['Incorrect number of arguments (supplied: ' num2str(nargin) '):\n Usage:       function [strand_collection, bundle_colours] = display_strands([dirname or strand_struct], [tube_divs], [bundle_colours],[bundle or strands to_plot],["bundle" or "strand" (which to plot)])\n']);
	end

	isotropic_regions = [];
	
	if (ischar(varargin{1}))
		[strand_collection, isotropic_regions] = load_strand_collection(varargin{1});
		collection_name = varargin{1};
	elseif (iscell(varargin{1}))

 		strand_collection = varargin{1};
    
    if size(strand_collection,1) == 0
      error('No strands in loaded collection');
    end
    
		collection_name = strand_collection{end,1};
    strand_collection(end,:) = [];
    
	else
		error('Incorrect data type for argument 1\n');    
	end   

	if (nargin >= 2) 
		bundle_colours = varargin{2};
	else
		bundle_colours = [];
	end

	if (nargin >= 3) 
		plot_indices = varargin{3};
	else
		plot_indices = [];
    end	
	
    if (nargin >= 4) 
		plot_indices_type = varargin{4};
	else
		plot_indices_type = 'bundle';
	end	
    
	if (nargin >= 5 && strcmp('line', varargin{5}))
		plot_3D = 0;
	else
		plot_3D = 1;
	end
	
	
	if (nargin == 6) 
		tube_divs = varargin{6};
	else
		
		if (size(strand_collection,1) <= 100)
			tube_divs = 32;
		elseif (size(strand_collection,1) <= 250)
			tube_divs = 16;	
		else		
			tube_divs = 6;
		end
	end 
	
	display_key = (size(bundle_colours,1) == 0);
		
	bundle_indices = [];

	if (size(bundle_colours,1) == 0)

		for strand_i = 1:size(strand_collection,1)

			bundle_indices = [bundle_indices; strand_collection{strand_i, 4}]; %#ok<AGROW>

		end

		bundle_indices = unique(bundle_indices);
    bundle_indices = sort(bundle_indices);

		bundle_colours = rand([(max(bundle_indices)+1), 3]);

	end
            

	if (display_key)
		display_colour_key(bundle_colours, bundle_indices);
	end
	
	
%	h = sfigure(gcf);
%	main_fig = sfigure(gcf);

	main_fig = figure();

	set(main_fig,'Units','normalized') 

	set(main_fig, 'Position', [0.3 0.25 0.4 0.5]);
	set(main_fig, 'DoubleBuffer', 'on');
	set(main_fig, 'Name', collection_name);
	
	cameratoolbar('Show');
	cameratoolbar('SetMode','orbit');

	clf;

  whitebg(main_fig,'black');
	hold on;

% 	new_bundle_i = [];

	tube_corners = tube_divs + 1;

	strand_count = 0;
	
	for (isotropic_region_i = 1:size(isotropic_regions,1))
		
		centre = isotropic_regions(isotropic_region_i, 1:3)';
		radius = isotropic_regions(isotropic_region_i, 4);
		
		num_faces = round(radius * 100);
		
		[X, Y, Z] = sphere(num_faces);
	
		X = X * radius + centre(1);
		Y = Y * radius + centre(2);
		Z = Z * radius + centre(3);
		
		h = surf(X, Y, Z);
		
		set(h,'facecolor', [0.5 0.5 0.5]);
		set(h, 'FaceAlpha', 0.5);
	end	

	for strand_i = 1:size(strand_collection,1)
        
%        dummy1 = strcmp(plot_indices_type,'bundle');
%        dummy2 = any(strand_collection{strand_i,4} == plot_indices);
%        dummy3 = strcmp(plot_indices_type,'strand');
%        dummy4 = any(strand_collection{strand_i,2} == plot_indices);
%        dummy5 = size(plot_indices,1) == 0;
    if (size(strand_collection{strand_i,1},1) < 2)
			disp(['Strand ' num2str(strand_i) ' omitted due to it only having ' num2str(size(strand_collection{strand_i,1},1)) ' control points.']);
			
		elseif ( (strcmp(plot_indices_type,'bundle') && any(strand_collection{strand_i,4} == plot_indices)) || (strcmp(plot_indices_type,'strand') && any(strand_collection{strand_i,2} == plot_indices)) || size(plot_indices,1) == 0) 

			strand = strand_collection{strand_i,1};
			strand_colour = bundle_colours(strand_collection{strand_i, 4}+1,:);
			
			
			if (plot_3D)
				
				save_strand_index = strand_collection{strand_i,2};
				strand_r = strand_collection{strand_i, 3}; 
				save_bundle_index = strand_collection{strand_i,4};


				num_control_points = size(strand,1);

				segs = strand(2:num_control_points,:) - strand(1:num_control_points-1,:);

				norm_segs = zeros(size(segs));

								
				length_segs = sqrt(dot(segs,segs,2));

				zero_length_segs = find(length_segs == 0.0);
				if (size(zero_length_segs,1) >= 1) 
					length_segs(zero_length_segs) = 0.0001;
					segs(zero_length_segs) = 0.0001;
					strand(zero_length_segs) = strand(zero_length_segs)  + 0.0001; 
				end

				norm_segs(:,1) = segs(:,1) ./ length_segs;
				norm_segs(:,2) = segs(:,2) ./ length_segs;
				norm_segs(:,3) = segs(:,3) ./ length_segs;

				consec_segs_align = dot(norm_segs(1:num_control_points-2,:), norm_segs(2:num_control_points-1,:),2);
						
				
			
				count = 0;
				while (1)
				 
					ref_vector = rand(1,3);					
				
					[t,n,b] = frame( strand(:,1),strand(:,2),strand(:,3), ref_vector);

					if (all(~isnan(n)))
						break;
					end
					
					if (count > 100) 
						error(['likely nan values in strand ' num2str(strand_i-1) '!']);
					end
					
					count = count + 1;
				end
					


				X = zeros(num_control_points, tube_corners);
				Y = zeros(num_control_points, tube_corners);
				Z = zeros(num_control_points, tube_corners);

				theta = 0:(2*pi/(tube_corners-1)):(2*pi);

				count = 0;

				for point_i = 1:num_control_points

					if (point_i == 1) 
						w = strand(1, :) + n(1,:);
						n_prime = n(1,:);
						b_prime = b(1,:);

					else




						mu = dot(t(point_i,:), strand(point_i,:) - w, 2) / dot(t(point_i,:), segs(point_i-1,:),2);

						w_proj = w + mu * segs(point_i-1, :);

						n_prime = w_proj - strand(point_i,:);

						n_prime = n_prime ./ norm(n_prime);

						b_prime = cross( t(point_i,:), n_prime);

						b_prime = b_prime ./ norm(b_prime);

						w = strand(point_i,:) + n_prime;

					end

					X(point_i,:) = strand(point_i, 1) + strand_r * ( n_prime(1,1) * cos(theta) + b_prime(1,1) * sin(theta));
					Y(point_i,:) = strand(point_i, 2) + strand_r * ( n_prime(1,2) * cos(theta) + b_prime(1,2) * sin(theta));
					Z(point_i,:) = strand(point_i, 3) + strand_r * ( n_prime(1,3) * cos(theta) + b_prime(1,3) * sin(theta));


				end

				h = surf(X,Y,Z);
				
% 
% 				if (bundle_i ~= (strand_collection{strand_i, 4}+1)) 
% 					bundle_i = (strand_collection{strand_i, 4}+1);
% 					strand_colour = bundle_colours(bundle_i,:);
% 					new_bundle_i = [new_bundle_i; bundle_i, strand_i];
% 
% 
% 				end

				set(h,'facecolor', strand_colour);
			else
				h = plot3(strand(:,1), strand(:,2), strand(:,3));
				set(h, 'Color', strand_colour);
			end

			strand_count = strand_count + 1;
			if (mod(strand_count,25) == 0)
				fprintf('.');
			end
    end


  end
  
  strand_collection{end+1,1} = collection_name;
  
	hold off;

	%set(gca, 'color', [0 0 0]);
	set(gcf, 'color', [0 0 0]);
% 	xlabel('X-axis', 'color', [1 1 1]);
% 	ylabel('Y-axis', 'color', [1 1 1]);
% 	zlabel('Z-axis', 'color', [1 1 1]);
	
	label_handle = get(gca,'xlabel');
	set(label_handle, 'string', 'X-axis', 'color', [1 1 1]);
% 	label_position = get(label_handle, 'position');
% 	set(label_handle, 'position', (label_position + [0,  0, 0]));
%   
	label_handle = get(gca,'ylabel');
	set(label_handle, 'string', 'Y-axis', 'color', [1 1 1]);
% 	label_position = get(label_handle, 'position');
% 	set(label_handle, 'position', (label_position + [0,  0, 0]));
   
	label_handle = get(gca,'zlabel');
	set(label_handle, 'string', 'Z-axis', 'color', [1 1 1]);
% 	label_position = get(label_handle, 'position');
% 	set(label_handle, 'position', (label_position + [0,  0, 0]));
   
	
	
	
	h = get (gca, 'children');
	daspect ([ 1 1 1 ]);
	if (plot_3D) 
		set(h,'edgecolor', 'none')
	end
	light           
	lighting gouraud;
	

%	figure(main_fig);

	fprintf('\n');
	disp(['Displaying ' num2str(strand_count) ' strands.']);
	
end


function [strand_collection, isotropic_regions] = load_strand_collection(dirname)

%  function strand_collection = load_strand_collection(dirname)
% 
%   load_strand_collection.m
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


	num_strands = 0;

	isotropic_regions = [];
	
	files = dir(dirname)'; %'

	if (size(files) == [1,0]) 
		error(['Could not load any strands from directory ' dirname ]);
	end

	for file = files

		if (~file.isdir)

			delimeters = [strfind(file.name, '_') strfind(file.name, '-') strfind(file.name, '.txt' )];

			if (length(delimeters) == 4 && strmatch('strand', file.name))

				num_strands = num_strands + 1;


				strand_collection{num_strands, 1} = load([dirname filesep file.name]);
				strand_collection{num_strands, 2} = str2num(file.name(delimeters(1) + 1 :delimeters(2) -1 ));
				strand_collection{num_strands, 3} = str2num(file.name(delimeters(3) + 2 :delimeters(4) -1 ));
				strand_collection{num_strands, 4} = str2num(file.name(delimeters(2) + 1 :delimeters(3) -1 ));


			end
			
			if (strcmp(file.name, 'isotropic_regions.txt'))
				isotropic_regions = load([dirname filesep file.name]);
				
			end
				
		end

	end
end



function [t,n,b]=frame(x,y,z,vec)

% FRAME Calculate a Frenet-like frame for a polygonal space curve
% [t,n,b]=frame(x,y,z,v) returns the tangent unit vector, a normal
% and a binormal of the space curve x,y,z. The curve may be a row or
% column vector, the frame vectors are each row vectors. 
%
% This function calculates the normal by taking the cross product
% of the tangent with the vector v; if v is chosen so that it is
% always far from t the frame will not twist unduly.
% 
% If two points coincide, the previous tangent and normal will be used.
%
% Written by Anders Sandberg, asa@nada.kth.se, 2005


	N=size(x,1);
	if (N==1)
	  x=x'; %'
	  y=y'; %'
	  z=z'; %'
	  N=size(x,1);
	end

	t=zeros(N,3);
	b=zeros(N,3);
	n=zeros(N,3);

	p=[x y z];

	for i=2:(N-1)
	  t(i,:)=(p(i+1,:)-p(i-1,:));
	  tl=norm(t(i,:));
	  if (tl>0)
		t(i,:)=t(i,:)/tl;
	  else
		t(i,:)=t(i-1,:);
	  end
	end

	t(1,:)=p(2,:)-p(1,:);
	t(1,:)=t(1,:)/norm(t(1,:));

	t(N,:)=p(N,:)-p(N-1,:);
	t(N,:)=t(N,:)/norm(t(N,:));

	for i=2:(N-1)
	  n(i,:)=cross(t(i,:),vec);
	  nl=norm(n(i,:));
	  if (nl>0)
		n(i,:)=n(i,:)/nl;
	  else
		n(i,:)=n(i-1,:);
	  end
	end

	n(1,:)=cross(t(1,:),vec);
	n(1,:)=n(1,:)/norm(n(1,:));

	n(N,:)=cross(t(N,:),vec);
	n(N,:)=n(N,:)/norm(n(N,:));

	for i=1:N
	  b(i,:)=cross(t(i,:),n(i,:));
	  b(i,:)=b(i,:)/norm(b(i,:));
	end
end
	
	
function h = sfigure(h)
% SFIGURE  Create figure window (minus annoying focus-theft).
%
% Usage is identical to figure.
%
% Daniel Eaton, 2005
%
% See also figure

	if nargin>=1 
		if ishandle(h)
			set(0, 'CurrentFigure', h);
		else
			h = figure(h);
		end
	else
		h = figure;
	end
end