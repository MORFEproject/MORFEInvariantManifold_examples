function out = curve
%
% Curve file of circle
%
    out{1}  = @curve_func;
    out{2}  = @defaultprocessor;
    out{3}  = @options;
    out{4}  = [];%@jacobian;
    out{5}  = [];%@hessians;
    out{6}  = [];%@testf;
    out{7}  = [];@userf;
    out{8}  = [];@process;
    out{9}  = [];@singmat;
    out{10} = [];@locate;
    out{11} = [];@init;
    out{12} = [];@done;
    out{13} = @adapt;
%---------------------------------------------------------  
function func = curve_func(arg);
  x = arg;
  func = x(1)^2+x(2)^2-1;
%---------------------------------------------------------
function option= options
  option = contset;
  varargout{1} = option;
%----------------------------------------------------------
function varargout = defaultprocessor(varargin)
  if nargin > 2
    s = varargin{3};
    varargout{3} = s;
  end
  
  % no special data
  varargout{2} = NaN;

  % all done succesfully
  varargout{1} = 0;
%----------------------------------------------------------  
function [res,x,v] = adapt(x,v)
res = []; % no re-evaluations needed
  
