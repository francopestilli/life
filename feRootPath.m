function p = feRootPath
% Returns the path to the Microtrack code.
%
% Example:
%   feP = feRootPath;
%
% See also: v_lifeExample
%
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.
p = which('feRootPath');

p = fileparts(p);

end