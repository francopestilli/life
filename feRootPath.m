function p = feRootPath
% Returns the path to the Microtrack code.
%
% Example:
%   feP = feRootPath;
%
% See also: v_lifeExample
%
% (c) Stanford VISTA Team

p = which('feRootPath');

p = fileparts(p);

end