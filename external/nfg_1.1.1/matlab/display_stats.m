function display_stats(dirname)

%	function display_stats(dirname)
%
%	args: 
%            dirname - name of the input directory (the output directory of a 'mri_sim' call with 'save_ext_stats' set to 1)
%	
%	
%
%	NB: To generate the extended statistics, set the input parameter 'save_ext_stats' to 1 and run 'mri_sim' on a strand collection. 
%         
%
%	display_stats
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



filename = [dirname '/stats/stats.txt'];

file = fopen(filename);

if (file == -1) 
	error([ 'Could not open file ' filename '!' ]);
end

line = fgetl(file);

disp(' ');
disp('---------');
disp(' ');

while (ischar(line)) 
    
    disp([line ' ']);
    line = fgetl(file);
end

disp('---------');
disp(' ');
