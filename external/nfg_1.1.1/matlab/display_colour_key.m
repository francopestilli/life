function display_colour_key(colours, bundle_indices)
%function display_colour_key(colours, bundle_indices)
%
%
%   Usage:	arg[1] = The colour specifications returned by 'display_strands'
%			arg[2] = The bundle colours to display
%
%   Takes a Mx3 colourspec matrix and displays it as an image to act as a
%   defacto legend for 'display_strands' 
%
%	display_colour_key
%	Numerical fibre Generator
%
%
%	Created by Tom Close on 17/07/07 
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

global colours_of_strands

if ~exist('colours') && ~isempty(colours_of_strands)
  
  colours = colours_of_strands;
  bundle_indices = [0:1:(size(colours,1)-1)]';
  
end

h = figure();
set(h,'Units','normalized');
set(h, 'Position', [0.0 0.05 0.05 0.9]);
colours2 = reshape(colours, size(colours,1),1,3);
% image(1, [0:1:(size(colours,1)-1)], colours2);
image(1, [0:1:(size(bundle_indices,1)-1)], colours2((bundle_indices+1),:,:));
set(gca, 'YTick', [0:1:(size(bundle_indices,1)-1)]);
set(gca, 'YTickLabel', bundle_indices);
set(gca, 'XTick', []);
set(get(gca,'YLabel'),'String','Bundle Index');