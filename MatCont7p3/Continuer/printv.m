function printv(varargin)

if nargin==1
  printconsole('( ');
  printconsole('%f ', varargin{1});
  printconsole(')\n');
else
  printconsole('%s = ( ', varargin{1});
  printconsole('%f ', varargin{2});
  printconsole(')\n');
end


%SD:prints a vector
