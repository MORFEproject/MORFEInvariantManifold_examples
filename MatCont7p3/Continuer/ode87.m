function varargout= ode87(ode,tspan,x0,options,varargin)
% ODE87  is a realization of explicit Runge-Kutta method. 
% Integrates a system of ordinary differential equations using
% 8-7 th order Dorman and Prince formulas.  See P.J. Prince & J.R. Dorman (1981) 
% High order embedded Runge-Kutta formulae. J.Comp. Appl. Math., Vol. 7. p.67-75.
%
% This is a 8th-order accurate integrator therefore the local error normally
% expected is O(h^9).  
% This requires 13 function evaluations per integration step.
%
% Some information about method can be found in
% Hairer, Norsett and Wanner (1993): Solving Ordinary Differential Equations. Nonstiff Problems. 
% 2nd edition. Springer Series in Comput. Math., vol. 8. 
%
%  Interface to program based on standart MATLAB ode-suite interface but
%  with some restriction. 
%   [T,Y] = ODE87(ODEFUN,TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates the
%   system of differential equations y' = f(t,y) from time T0 to TFINAL with
%   initial conditions Y0. Function ODEFUN(T,Y) must return a column vector
%   corresponding to f(t,y). Each row in the solution array Y corresponds to
%   a time returned in the column vector T. 
%   
%   [T,Y] = ODE87(ODEFUN,TSPAN,Y0,OPTIONS) solves as above with default
%   integration properties replaced by values in OPTIONS, an argument created
%   with the ODESET function. See ODESET for details. Commonly used options 
%   are scalar relative error tolerance 'RelTol' (1e-6 by default).
%   
%   [T,Y] = ODE87(ODEFUN,TSPAN,Y0,OPTIONS,P1,P2...) passes the additional
%   parameters P1,P2,... to the ODE function as ODEFUN(T,Y,P1,P2...), and to
%   all functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
%   no options are set.   
%
%   Example    
%         [t,y]=ode87(@vdp1,[0 20],[2 0]);   
%         plot(t,y(:,1));
%     solves the system y' = vdp1(t,y), using the default relative error
%     tolerance 1e-3 and the default absolute tolerance of 1e-6 for each
%     component, and plots the first component of the solution. 
%
% --------------------------------------------------------------------
% Copyright (C) 2003, Govorukhin V.N.
% This file is intended for use with MATLAB and was produced for MATDS-program
% http://www.math.rsu.ru/mexmat/kvm/matds/
% ODE87 is free software. ODE87 is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY. 



% The coefficients of method
c_i=  [ 1/18, 1/12, 1/8, 5/16, 3/8, 59/400, 93/200, 5490023248/9719169821, 13/20, 1201146811/1299019798, 1, 1]';

a_i_j = [ 1/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
          1/48, 1/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
          1/32, 0, 3/32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
          5/16, 0, -75/64, 75/64, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
          3/80, 0, 0, 3/16, 3/20, 0, 0, 0, 0, 0, 0, 0, 0; 
          29443841/614563906, 0, 0, 77736538/692538347, -28693883/1125000000, 23124283/1800000000, 0, 0, 0, 0, 0, 0, 0;
          16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777, 545815736/2771057229, -180193667/1043307555, 0, 0, 0, 0, 0, 0;
          39632708/573591083, 0, 0, -433636366/683701615, -421739975/2616292301, 100302831/723423059, 790204164/839813087, 800635310/3783071287, 0, 0, 0, 0, 0;
          246121993/1340847787, 0, 0, -37695042795/15268766246, -309121744/1061227803, -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 123872331/1001029789, 0, 0, 0, 0;
         -1028468189/846180014, 0, 0, 8478235783/508512852, 1311729495/1432422823, -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649, -45442868181/3398467696, 3065993473/597172653, 0, 0, 0;
          185892177/718116043, 0, 0, -3185094517/667107341, -477755414/1098053517, -703635378/230739211, 5731566787/1027545527, 5232866602/850066563, -4093664535/808688257, 3962137247/1805957418, 65686358/487910083, 0, 0;
          403863854/491063109, 0, 0, -5068492393/434740067, -411421997/543043805, 652783627/914296604, 11173962825/925320556, -13158990841/6184727034, 3936647629/1978049680, -160528059/685178525, 248638103/1413531060, 0, 0]';

 b_8 = [ 14005451/335480064, 0, 0, 0, 0, -59238493/1068277825, 181606767/758867731,   561292985/797845732,   -1041891430/1371343529,  760417239/1151165299, 118820643/751138087, -528747749/2220607170,  1/4]';

 b_7 = [ 13451932/455176623, 0, 0, 0, 0, -808719846/976000145, 1757004468/5645159321, 656045339/265891186,   -3867574721/1518517206,   465885868/322736535,  53011238/667516719,                  2/45,    0]';


pow = 1/8; % power for step control

% Check inputs
if nargin < 5
   varargin={};
end;
if nargin < 4
  options = [];
  if nargin < 3
     error('Not enough input arguments.  See ODE87.');
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
end;
% initial step size
FcnHandlesUsed = isa(ode,'function_handle');
soloutRequested = (FcnHandlesUsed & (nargout==1));

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, args, ...
          options, atol, rtol, threshold, normcontrol, normy, hmax, htry, htspan]  ...
         = odearguments(FcnHandlesUsed,'ode87', ode, tspan, x0, options,  ...
                        soloutRequested, varargin);
one2neq = (1:neq)';
h=odeget(options,'InitialStep');
if isempty(h)
   h = (tspan(2) - tspan(1))/50;
   if h>0.1
      h=0.1;
   end
   if h>hmax 
      h = hmax;
   end;
end;

% Output ODEction checking and output parameters
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

%Initialization
t0 = tspan(1);
tfinal = tspan(2);
t = t0;

% Minimal step size
hmin = 16*eps*abs(t);


% constant for step rejection
reject = 0;

x = x0(:);          % start point
f = x*zeros(1,13);  % array f for RHS calculation
tau = tol * max(norm(x,'inf'), 1);  % accuracy

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
      if (t + h) > tfinal 
         h = tfinal - t; 
      end;

      nsteps =nsteps + 1;

% Compute the RHS for step of method
      f(:,1) = feval(ode,t,x,args{:});
      for j = 1: 12,
          f(:,j+1) = feval(ode, t+c_i(j)*h, x+h*f*a_i_j(:,j),args{:});
      end;
      nfevals= nfevals+13;

% Two solution 
      sol2=x+h*f*b_8;
      sol1=x+h*f*b_7;

% Truncation error 
      error_1 = norm(sol1-sol2);

% Estimate the error and the acceptable error

      Error_step = norm(error_1,'inf');
      tau = tol*max(norm(x,'inf'),1.0);

% Update the solution only if the error is acceptable

      if Error_step <= tau
         tnew = t + h;
         xnew = sol2; 
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
      end;

% Step control
      if Error_step == 0.0
         Error_step = eps*10.0;
      end;
      h = min(hmax, 0.9*h*(tau/Error_step)^pow);
      if (abs(h) <= eps)|isnan(Error_step) 
         if reject == 0
            disp('Warning!!! ode87. Step is very small!!!');
            h = eps * 100;
            return
         else
            disp('Error in ode87. Step is too small.');
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

  end;

   if (t < tfinal)
    varargout{1} = tout(1:nout);
    varargout{2} = yout(1:nout,:);
    if haveoutfun
        feval(outfun,[],[],'done',outputArgs1{:});
    end
      disp('Error in ODE87...')
   end;

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