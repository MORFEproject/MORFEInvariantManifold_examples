function [neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, args, ...
          options, atol, rtol, threshold, normcontrol, normy, hmax, htry, htspan]    ...
         = odearguments(FcnHandlesUsed, solver, ode, tspan, y0, options, ...
                        soloutRequested, extras)
%ODEARGUMENTS  Helper function that processes arguments for all ODE solvers.
%
%   See also ODE113, ODE15S, ODE23, ODE23S, ODE23T, ODE23TB, ODE45.

%   Mike Karr, Jacek Kierzenka
%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.12 $  $Date: 2002/04/08 20:26:55 $

if FcnHandlesUsed  % function handles used
  msg = ['When the first argument to ', solver,' is a function handle, '];
  if isempty(tspan) | isempty(y0) 
    error([msg 'the tspan and y0 arguments must be supplied.']);
  end      
  if length(tspan) < 2
    error([msg 'the tspan argument must have at least two elements.']);
  end  
  htspan = abs(tspan(2) - tspan(1));
  if soloutRequested  % output sol structure 
    tspan = [tspan(1); tspan(end)];
  else
    tspan = tspan(:);
  end  
  ntspan = length(tspan);
  t0 = tspan(1);  
  next = 2;       % next entry in tspan
  tfinal = tspan(end);     
  args = extras;                 % use f(t,y,p1,p2...) 

else  % ode-file used
  if ~( isempty(options) | isa(options,'struct') )
    if (length(tspan) == 1) & (length(y0) == 1) & (min(size(options)) == 1)
      tspan = [tspan; y0];
      y0 = options;
      options = [];
      warning('MATLAB:odearguments:ObsoleteSyntax',['Obsolete syntax.  Use ' ...
               '%s(fun,tspan,y0,...) instead.'],solver);
    else
      error(['Correct syntax is ', solver, '(', funstring(ode), ...
             ',tspan,y0,options).']);      
    end
  end   
  % Get default tspan and y0 from the function if none are specified.
  if isempty(tspan) | isempty(y0) 
    if exist(ode)==2 & ( nargout(ode)<3 & nargout(ode)~=-1 ) 
      % nargout M-files only
      msg = sprintf('Use %s(%s,tspan,y0,...) instead.',solver,funstring(ode));
      error(['No default parameters in ' funstring(ode) '.  ' msg]);
    end
    [def_tspan,def_y0,def_options] = feval(ode,[],[],'init',extras{:});
    if isempty(tspan)
      tspan = def_tspan;
    end
    if isempty(y0)
      y0 = def_y0;
    end
    options = odeset(def_options,options);
  end  
  tspan = tspan(:);
  ntspan = length(tspan);
  if ntspan == 1    % Integrate from 0 to tspan   
    t0 = 0;          
    next = 1;       % Next entry in tspan.
  else              
    t0 = tspan(1);  
    next = 2;       % next entry in tspan
  end
  htspan = abs(tspan(next) - t0);
  tfinal = tspan(end);   
  
  % The input arguments of f determine the args to use to evaluate f.
  if (exist(ode)==2)                 % M-file
    if (nargin(ode) == 2)           
      args = {};                   % f(t,y)
    else
      args = [{''} extras];        % f(t,y,'',p1,p2...)
    end
  else  % MEX-files, etc.
    try 
      args = [{''} extras];        % try f(t,y,'',p1,p2...)     
      dummy = feval(ode,tspan(1),y0(:),args{:});   
    catch
      lasterr('');
      args = {};                   % use f(t,y) only
    end
  end
end

y0 = y0(:);
neq = length(y0);

% Test that tspan is internally consistent.
if t0 == tfinal
  error('The last entry in tspan must be different from the first entry.');
end
tdir = sign(tfinal - t0);
if any( tdir*diff(tspan) <= 0 )
  error('The entries in tspan must strictly increase or decrease.');
end

f0 = feval(ode,t0,y0,args{:});
[m,n] = size(f0);

if n > 1
  error([funstring(ode) ' must return a column vector.'])
elseif m ~= neq
  msg = sprintf('an initial condition vector of length %d.',m);
  error(['Solving ' funstring(ode) ' requires ' msg]);
end

% Get the error control options, and set defaults.
rtol = odeget(options,'RelTol',1e-3,'fast');
if (length(rtol) ~= 1) | (rtol <= 0)
  error('RelTol must be a positive scalar.');
end
if rtol < 100 * eps 
  rtol = 100 * eps;
  warning('MATLAB:odearguments:RelTolIncrease', ...
          'RelTol has been increased to %g.',rtol)
end
atol = odeget(options,'AbsTol',1e-6,'fast');
if any(atol <= 0)
  error('AbsTol must be positive.');
end
normcontrol = strcmp(odeget(options,'NormControl','off','fast'),'on');
if normcontrol
  if length(atol) ~= 1
    error('Solving with NormControl ''on'' requires a scalar AbsTol.');
  end
  normy = norm(y0);
else
  if (length(atol) ~= 1) & (length(atol) ~= neq)
    error(sprintf(['Solving %s requires a scalar AbsTol, or a vector' ...
                   ' AbsTol of length %d'],funstring(ode),neq)); 
  end
  atol = atol(:);
  normy = [];
end
threshold = atol / rtol;

% By default, hmax is 1/10 of the interval.
hmax = min(abs(tfinal-t0), abs(odeget(options,'MaxStep',0.1*(tfinal-t0),'fast')));
if hmax <= 0
  error('Option ''MaxStep'' must be greater than zero.');
end
htry = abs(odeget(options,'InitialStep',[],'fast'));
if ~isempty(htry) & (htry <= 0)
  error('Option ''InitialStep'' must be greater than zero.');
end
