function long_fibers = feFindLongFibers(fibers,min_length)
%
% Find indices to fibers longer than a minimum length.
% 
% long_fibers = feFindLongFibers(fg,min_length)
%
% INPUTS: 
%        fibers    - A cell array of fibers. e.g., fg.fibers.
%        min_lenth - Minimum length in nodes (number of x,y,z
%                    coordinates).
% 
% OUTPUTS:
%        long_fibers - A vector of zeros and ones the same length of
%                      fibers. The vector contains a 1 for each long fiber,
%                      0 otherwise.
%
% SEE ALSO: feConnectomePreprocess.m
%
% Franco (c)  2012 Stanford Vista Team.

% Fiber length is defined in nodes (e.g., number of xyz coordinates in a
% fiber).
if ~exist('min_length','var') || isempty(min_length)
  min_length = 10; % Default minimum length is 10 nodes
  % fprintf('[%s] Discarding fibers with fewer than %i nodes',mfilename, min_length)
end

% Find fibers that are longer than min Length
long_fibers = cellfun(@length, fibers) >= min_length;

end % End main function