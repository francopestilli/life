function p = feRootPath
% Returns the path to the Microtrack code.
%
% Example:
%   feP = feRootPath;
%
% See also: v_lifeExample
%
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com.
p = which('feRootPath');

p = fileparts(p);

end
