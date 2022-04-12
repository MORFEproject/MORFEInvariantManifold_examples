function debug(varargin)
%
% Debug routine, initially for testing purposes.
% It also provides some extra output in the Matlab window during runs.
%

if nargin == 0
  error('debug needs at least one parameter');
end

if nargin == 1
  fprintf(varargin{1});
else
  fprintf(varargin{1}, varargin{2:nargin});
end
