function val = contget(options, name, default)
%
% val = CONTGET(options, name, default)
%
% Returns value of an option
%
% Input : options: option vector created with CONTSET
%         name   : name of option
%         default: default value if option not set
% Output: val: value of option or default (input)

opt = contidx(name);

if isempty(opt)
  msg = sprintf('Unrecognized option ''%s'' ', name);
  error(msg);
else
  eval(['val = options.' opt ';']);
  if isempty(val)
    val = default;
  end
end

%SD:gets value of option
