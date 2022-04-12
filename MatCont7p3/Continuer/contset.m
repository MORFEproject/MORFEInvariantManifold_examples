function options = contset(options, name, value)
%
% options = CONTSET(options, optname, value)
%

if nargin == 0
  allopt = contidx;
  [m,n] = size(allopt);
  for i=1:m
    eval(['options.' allopt(i,:) '= [];']);
  end
  return;
end

if ~isstr(name)
  msg = sprintf('2nd argument must be a string');
  error(msg);
end

opt = contidx(name);

if isempty(opt)
  msg = sprintf('Unrecognized continuation option ''%s'' ', name);
  error(msg);
else
  eval(['options.' opt '= value;']);
end


%SD:sets value for option
