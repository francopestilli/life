function display_subvoxels(varargin)

%	function display_subvoxels(subvoxels, [slice_dim], [slice_num])
%
%	args: 
%            dirname - name of the input directory (the output directory of a 'mri_sim' call with 'save_subvoxels' set to 1)
%            slice_dim - the out of plane slice dimension.
%			 slice_num - the number of the slice to display.  If left blank
%						 all slices will be displayed in sequence (press
%						 any key to show the next slice)
%	
%	displays: A slice by slice representation of the plotted subvoxel orientations.  
%
%
%	NB: To generate the appropriate image, set the input parameter 'save_subvoxels' to 1 and run 'mri_sim' on a strand collection. 
%         
%
%	display_subvoxels
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

	if (nargin < 1 || nargin > 3) 
	   error('Incorrect number of arguments:\n Usage:       function display_subvoxels(subvoxels, [slice_dim], [slice_num])\n');
	end

	if (class(varargin{1}) == 'char')
		subvoxels = load_subvoxels(varargin{1});
	else
		error('Incorrect data type for argument 1\n');    
	end

	if (nargin >= 2) 
		slice_dim = varargin{2};
	else
		slice_dim = 1;
	end        

	if (nargin >= 3) 
		slice_num = varargin{3};
	else
		slice_num = [];
	end

	h = figure;

	set(h,'Units','normalized') 
	set(h, 'Position', [0.2 0.2 0.5 0.8]);

	subvoxels(find(subvoxels == -99)) = 0.75;
	
	
	
	if (slice_num ~= 0)

		
		if (slice_dim == 1)
			slice = squeeze(abs(subvoxels(slice_num,:,:,:)));
		elseif (slice_dim == 2)
			slice = squeeze(abs(subvoxels(:,slice_num,:,:)));
		elseif (slice_dim == 3)
			slice = squeeze(abs(subvoxels(:,:,slice_num,:)));
		end

		slice = permute(slice, [2 1 3]);	
		slice = slice(end:-1:1,end:-1:1,:);
		
		image(slice);
		daspect([1 1 1]);
		
		title(['Slice ' num2str(slice_num) ' of dim ' num2str(slice_dim)]);
	else


		for (i=[1:size(subvoxels,slice_dim)])

			if (slice_dim == 1)
				slice = squeeze(abs(subvoxels(i,:,:,:)));
			elseif (slice_dim == 2)
				slice = squeeze(abs(subvoxels(:,i,:,:)));
			elseif (slice_dim == 3)
				slice = squeeze(abs(subvoxels(:,:,i,:)));
			end

			slice = permute(slice, [2 1 3]);	
			slice = slice(end:-1:1,end:-1:1,:);
			
			image(slice);			
			daspect([1 1 1]);			
			title(['Slice ' num2str(i) ' of dim ' num2str(slice_dim)]);
			pause;

		end

		close(h);
	end	
	
	

end

function subvoxels = load_subvoxels(dirname)
% function subvoxels = load_subvoxels(filename)
%	args: 
%            dirname - name of the input directory (the output directory of a 'mri_sim' call with 'save_subvoxels' set to 1)


	filename = [dirname '/subvoxels/subvoxels.dat'];

	file = fopen(filename);

	if (file == -1) 
		error([ 'Could not open file ' filename '!' ]);
	end

	elems = fread(file, 'double');

	num_elems = size(elems,1);

	grid_length = nthroot( num_elems / 3, 3 );

	if (grid_length ~= round(grid_length))
		error('The number elements of the orientation plot cannot be evenly split in 3 dimensions (i.e. do not have an integer cube root)'); 
	end


	disp(['Grid dimensions: ' num2str(grid_length) ', ' num2str(grid_length) ', ' num2str(grid_length)]);

	subvoxels = reshape(elems, [3, grid_length, grid_length, grid_length]);
	subvoxels = shiftdim(subvoxels, 1);
end

