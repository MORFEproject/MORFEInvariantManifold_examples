function varargout = ode78(ode,tspan,x0,options,varargin)
% ODE78  is a realization of explicit Runge-Kutta method. 
% Integrates a system of ordinary differential equations using
% 7 th order Fehlberg formulas.  See Fehlberg (1969) 
% Classical fifth-, sixth-, seventh-, and eighth order Runge-Kutta formulas
% with step size control. Computing, Vol.4, p.93-106
%
% This is a 7-8th-order accurate integrator therefore the local error normally
% expected is O(h^9).  
% This requires 13 function evaluations per integration step.
%
% Some information about method can be found in
% Hairer, Norsett and Wanner (1993): Solving Ordinary Differential Equations. Nonstiff Problems. 
% 2nd edition. Springer Series in Comput. Math., vol. 8. 
%
%  Interface to program based on standart MATLAB ode-suite interface but
%  with some restriction. 
%   [T,Y] = ODE78(ODEFUN,TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates the
%   system of differential equations y' = f(t,y) from time T0 to TFINAL with
%   initial conditions Y0. Function ODEFUN(T,Y) must return a column vector
%   corresponding to f(t,y). Each row in the solution array Y corresponds to
%   a time returned in the column vector T. 
%   
%   [T,Y] = ODE78(ODEFUN,TSPAN,Y0,OPTIONS) solves as above with default
%   integration properties replaced by values in OPTIONS, an argument created
%   with the ODESET function. See ODESET for details. Commonly used options 
%   are scalar relative error tolerance 'RelTol' (1e-6 by default).
%   
%   [T,Y] = ODE78(ODEFUN,TSPAN,Y0,OPTIONS,P1,P2...) passes the additional
%   parameters P1,P2,... to the ODE function as ODEFUN(T,Y,P1,P2...), and to
%   all functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
%   no options are set.   
%
%   Example    
%         [t,y]=ode78(@vdp1,[0 20],[2 0]);   
%         plot(t,y(:,1));
%     solves the system y' = vdp1(t,y), using the default relative error
%     tolerance 1e-3 and the default absolute tolerance of 1e-6 for each
%     component, and plots the first component of the solution. 
%
% --------------------------------------------------------------------
% Copyright (C) 2003, Govorukhin V.N.
% This file is intended for use with MATLAB and was produced for MATDS-program
% http://www.math.rsu.ru/mexmat/kvm/matds/
% ODE78 is free software. ODE78 is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY. 

% The Fehlberg coefficients:
c_i = [ 2./27., 1/9, 1/6, 5/12, 0.5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1 ]';
a_i_j  = [ 2/27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
          1/36, 1/12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
          1/24, 0, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
          5/12, 0, -25/16, 25/16, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
          0.05, 0, 0, 0.25, 0.2, 0, 0, 0, 0, 0, 0, 0, 0 ;
          -25/108, 0, 0, 125/108, -65/27, 125/54, 0, 0, 0, 0, 0, 0, 0 ;
          31/300, 0, 0, 0, 61/225, -2/9, 13/900, 0, 0, 0, 0, 0, 0 ;
          2, 0, 0, -53/6, 704/45, -107/9, 67/90, 3, 0, 0, 0, 0, 0 ;
          -91/108, 0, 0, 23/108, -976/135, 311/54, -19/60, 17/6, -1/12, 0, 0, 0, 0 ;
          2383/4100, 0, 0, -341/164, 4496/1025, -301/82, 2133/4100, 45/82, 45/164, 18/41, 0, 0, 0 ;
          3/205, 0, 0, 0, 0, -6/41, -3/205, -3/41, 3/41, 6/41, 0, 0, 0 ;
          -1777/4100, 0, 0, -341/164, 4496/1025, -289/82, 2193/4100, 51/82, 33/164, 12/41, 0, 1, 0 ]';
b_7  = [ 0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0, 41/840, 41/840]';
b_8  = [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1 ]';



pow = 1/8; % constant for step control


% Check inputs
if nargin < 5
   varargin={};
end;
if nargin < 4
  options = [];
  if nargin < 3
     error('Not enough input arguments.  See ODE78.');
  end
end;
%Stats
nsteps  = 0;
nfailed = 0;
nfevals = 0;

% Maximal step size
hmax=odeget(options,'MaxStep');
if isempty(hmax)
   hmax = (tspan(2) - tspan(1))/2.5;
end

% initial step size
FcnHandlesUsed = isa(ode,'function_handle');
soloutRequested = (FcnHandlesUsed & (nargout==1));

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, args, ...
          options, atol, rtol, threshold, normcontrol, normy, hmax, htry, htspan]  ...
         = odearguments(FcnHandlesUsed,'ode78', ode, tspan, x0, options,soloutRequested, varargin);

one2neq = (1:neq)';
h=odeget(options,'InitialStep');
if isempty(h)
   h = (tspan(2) - tspan(1))/50;
   if h>0.1
      h=0.1;
   end;
   if h>hmax 
      h = hmax;
   end;
end;

% Output function checking and output parameters
if nargout > 0
  outfun = odeget(options,'OutputFcn',[],'fast');
else
  outfun = odeget(options,'OutputFcn',@odeplot,'fast');
end
if isempty(outfun)
  haveoutfun = 0;
else
  haveoutfun = 1;
  outputs = odeget(options,'OutputSel',one2neq,'fast');
  if isa(outfun,'function_handle')
    outputArgs = [{''},varargin];
    outputArgs1 = varargin;
  else       % With v5 syntax do not pass additional input arguments to outfun  
    outputArgs = {};      
    outputArgs1 = {};
  end  
end


%  A relative error tolerance that applies to all components of the solution vector. 
tol=odeget(options,'RelTol');
refine = odeget(options,'Refine',1,'fast');

if isempty(tol)
   tol = 1.e-6;
end;

% Initialization
t0 = tspan(1);
tfinal = tspan(2);
t = t0;

% Minimal step size
hmin = 16*eps*abs(t);

% constant for step rejection
reject = 0;


x = x0(:);          % start point
f = x*zeros(1,13);  % array f for RHS calculation
tau = tol * max(norm(x,'inf'), 1);%accuracy

% Initial output
if haveoutfun
  feval(outfun,[t tfinal],x(outputs),'init',outputArgs1{:});
end
if nargout > 0
    % alloc in chunks
  chunk = min(max(100,50*refine),floor((2^13)/neq));
  tout = zeros(chunk,1);
  yout = zeros(chunk,neq);
  nout = 1;
  tout(nout) = t;
  yout(nout,:) = x.';
end

% The main loop

while (t < tfinal) & (h >= hmin)
    if t + h > tfinal
        h = tfinal - t;        
    end
    
    % Compute RHS for step of method 
    f(:,1) = feval(ode,t,x,args{:});
    for j = 1: 12,
        f(:,j+1) = feval(ode, t+c_i(j)*h, x+h*f*a_i_j(:,j),args{:});
    end;
    nfevals = nfevals + 13;

    
    % Error on step
    error_1 = h*41/840*f*b_8;
    
    % Estimate the error and the acceptable error
    Error_step = norm(error_1,'inf');
    tau = tol*max(norm(x,'inf'),1.0);

    % Update and output the solution only if the error is acceptable
    if Error_step <= tau
        %nfailed = 0;
        tnew = t + h;
        xnew = x + h*f*b_7;  % this integrator uses local extrapolation
        if nargout > 0
             oldnout = nout;
             % computed points, no refinement
             nout = nout + 1;
             if nout > length(tout)
                 tout = [tout; zeros(chunk,1)];
                 yout = [yout; zeros(chunk,neq)];
             end
             tout(nout) = tnew;
             yout(nout,:) = xnew.';
         end   
         if haveoutfun
             i = oldnout+1:nout;
             if ~isempty(i) & (feval(outfun,tout(i),yout(i,outputs).',outputArgs{:}) == 1)
                 feval(outfun,[],[],'done',outputArgs1{:});                 
                 varargout{1} = tout(1:nout);
                 varargout{2} = yout(1:nout,:);
                 return;
             end
         end
         reject = reject - 1;
         t=tnew;
         x=xnew;
     else 
         nfailed = nfailed + 1;
         reject = 1;
     end
     
     % Update the step size
     if Error_step == 0.0
         Error_step = eps*10.0;
     end
     h = min(hmax, 0.8*h*(tau/Error_step)^pow);
     if (abs(h) <= eps) |isnan(Error_step)
         if reject == 0
             disp('Warning!!! ode78. Step is very small!!!');
             h = eps * 100;
         else
             disp('Error in ode78. Step is too small.');
             if nargout > 0
                 varargout{1} = tout(1:nout);
                 varargout{2} = yout(1:nout,:);
                 varargout{end+1} = [nsteps;nfailed;nfevals];
             end
             if haveoutfun
                 feval(outfun,[],[],'done',outputArgs1{:});
             end
             return;
         end;
     end;
     nsteps = nsteps + 1;
 end;

 if (t < tfinal)
    varargout{1} = tout(1:nout);
    varargout{2} = yout(1:nout,:);
    if haveoutfun
        feval(outfun,[],[],'done',outputArgs1{:});
    end
    disp('Error in ODE78')
    return
 end

if haveoutfun
      feval(outfun,[],[],'done',outputArgs1{:});
end
if nargout > 0
    varargout{1} = tout(1:nout);
    varargout{2} = yout(1:nout,:);
    varargout{end+1} = [nsteps;nfailed;nfevals];
end
printstats = strcmp(odeget(options,'Stats','off','fast'),'on');
stats = struct('nsteps',nsteps,'nfailed',nfailed,'nfevals',nfevals);

if printstats
    fprintf('%g successful steps\n', stats.nsteps);
    fprintf('%g failed attempts\n', stats.nfailed);
    fprintf('%g function evaluations\n', stats.nfevals);
end


