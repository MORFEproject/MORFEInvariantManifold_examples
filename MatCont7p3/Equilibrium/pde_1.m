function out = pde_1

    out{1}  = @curve_func;
    out{2}  = @defaultprocessor;
    out{3}  = @options;
    out{4}  = @jacobian;
    out{5}  = @hessians;
    out{6}  = @testf;
    out{7}  = @userf;
    out{8}  = @process;
    out{9}  = @singmat;
    out{10} = @locate;
    out{11} = @init;
    out{12} = @done;
    out{13} = @adapt;
return
%---------------------------------------------------
function func = curve_func(arg)
global eds 
  [x,p] = rearr(arg); p = num2cell(p);
  func = feval(eds.func, 0, x, p{:});
%----------------------------------------------------  
function jac = jacobian(varargin)
global eds 
    xo = varargin{1}; [x,p] = rearr(xo); p = num2cell(p);
    jac = [cjac(eds.func,eds.Jacobian,x,p,eds.ActiveParams) cjacp(eds.func,eds.JacobianP,x,p,eds.ActiveParams)];   
%----------------------------------------------------
function hess = hessians(varargin)
global eds 
  x = varargin{1};
  [x,p] = rearr(x); p = num2cell(p);
  h = chess(eds.func,eds.Jacobian,eds.Hessians,x,p,eds.ActiveParams);
  hp=chessp(eds.func,eds.Jacobian,eds.HessiansP,x,p,eds.ActiveParams);
  hess(:,:,size(x1,1)) = 0;
  for i=1:size(x,1)
      hess(:,:,i) = [h(:,:,i) hp(:,i)];
  end
%--------------------------------------------------
function varargout= defaultprocessor(varargin)
global cds eds
  if nargin > 2
    s = varargin{3};
    varargout{3} = s;
  end

  % compute eigenvalues?
  if (cds.options.Eigenvalues==1)
      xo = varargin{1}; [x,p] = rearr(xo); p = num2cell(p);
      varargout{2} = eig(cjac(eds.func,eds.Jacobian,x,p,eds.ActiveParams));
  else
      varargout{2} = nan;
  end  
  % all done succesfully
  varargout{1} = 0;
%---------------------------------------------------  
function option = options
global eds
  option = contset;

  % Check for symbolic derivatives in odefile
  
  symord = 0; 
  if ~isempty(eds.Jacobian), symord = 1; end
  if ~isempty(eds.Hessians), symord = 2; end
  if ~isempty(eds.Der3), symord = 3; end
  if ~isempty(eds.Der4), symord = 4; end
  if ~isempty(eds.Der5), symord = 5; end

  option = contset(option, 'SymDerivative', symord);
  option = contset(option, 'Singularities',1);
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Adapt', 0);
  option = contset(option, 'Locators', [0 0]);

  symordp = 0;
  if ~isempty(eds.JacobianP), symordp = 1; end
  if ~isempty(eds.HessiansP),  symordp = 2; end
  option = contset(option, 'SymDerivativeP', symordp);

%------------------------------------------------------
function [out, failed] = testf(id, x, v)
global cds 

if any(ismember(id,[1 2]))
    J=cjac(cds.curve_func,cds.curve_jacobian,x,[]);
end

out(2) = 0;
failed = [];

for i=id
  lastwarn('');
  
  switch i
  case 1 % BP
    % Jacobian extended with bordering vectors v and w
    B = [J; v'];
    out(1) = det(B);          
  case 2 % LP
    out(2) = v(end); 
    
  otherwise
    error('No such testfunction');
  end  
  if ~isempty(lastwarn)
    failed = [failed i];
  end
end
%--------------------------------------------------------
function [failed,s] = process(id, x, v, s)
global cds 
ndim = cds.ndim;

% WM: Removed SL array
fprintf('label = %s, x = ', s.label); printv(x);

switch id
  case 1 % BP
    s.msg  = sprintf('Branch point');
  case 2 % LP
    s.data.a=a_lp(x);
    fprintf('a=%d\n',s.data.a);
    s.msg  = sprintf('Limit point');
end

% Compute eigenvalues for every singularity
J=cjac(cds.curve_func,cds.curve_jacobian,x,[]);
if ~issparse(J)
  [v,d]=eig(J(:,1:ndim-1));
else
  opt.disp=0;
  % WM: fixed a bug (incompatability between MatLab 6.0 and 5.5?)
  [v,d]=eigs(J(:,1:ndim-1),min(6,ndim-1),'lm',opt);
end

s.data.evec = v;
s.data.eval = diag(d)';

failed = 0;
%------------------------------------------------------------
function [S,L] = singmat    
% 0: testfunction must vanish
% 1: testfunction must not vanish
% everything else: ignore this testfunction

  S = [  0 8
         1 0 ];

  L = [ 'BP'; 'LP' ];
%-------------------------------------------------------------
function [x,v] = locate(id, x1, v1, x2, v2)
switch id
  case 1
    [x,v] = locateBP(id, x1, v1, x2, v2);
  otherwise
    error('No locator defined for singularity %d', id);
end
%--------------------------------------------------------------
function varargout = init(varargin)
  x = varargin{1};
  v = varargin{2};
  WorkspaceInit(x,v);

  % all done succesfully
  varargout{1} = 0;
%--------------------------------------------------------------
function done
  WorkspaceDone;
%--------------------------------------------------------------
function [res,x,v] = adapt(x,v)
res = []; % no re-evaluations needed

% ---------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------


function [x,p] = rearr(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x) and parameters (p)
global cds eds

nap = length(eds.ActiveParams);
ncoo = cds.ndim-nap;

p = eds.P0;
p(eds.ActiveParams) = x0((ncoo+1):end);
x = x0(1:ncoo);

% ---------------------------------------------------------------

function [x,v] = locateBP(id, x1, v1, x2, v2)
global cds

ndim = cds.ndim;

J = cjac(cds.curve_func,cds.curve_jacobian,x1,[]);
if ~issparse(J)
  [v,d]=eig(J(:,1:ndim-1));
else
  opt.disp=0;
  [v,d]=eigs(J(:,1:ndim-1), 'SM', opt);
end

p = v(:,find(min(abs(diag(d)))));
b = 0;
x = x1;
i = 0;

u = [x; b; p];

[A,f]=locjac(x,b,p);

while i < cds.options.MaxCorrIters
  
  du = A\f;
  u = u - du;

  x = u(1:ndim);
  b = u(ndim+1);
  p = u(ndim+2:2*ndim);

  [A,f]=locjac(x,b,p);

  % WM: VarTol and FunTol were switched
  if norm(du) < cds.options.VarTolerance & norm(f) < cds.options.FunTolerance break; end

  i = i+1;
end

v = 0.5*(v1+v2);

% ---------------------------------------------------------------

function [A, f] = locjac(x, b, p)
% A = jac of system
% f = system evaluated at (x,b,p)
global cds 

ndim = cds.ndim;

II = eye(ndim-1);
J = cjac(cds.curve_func,cds.curve_jacobian,x,[]);
H = chess(cds.curve_func,cds.curve_jacobian,cds.curve_hessians,x,[]);

F1 = [J, p, b*II];

for j=1:ndim
  for k=j:ndim
    F21(j,k) = H(:,j,k)'*p;
    F21(k,j) = F21(j,k);
  end
end

F22 = zeros(ndim,1);
F23 = J';

F3 = [zeros(1,ndim), 0, 2*p'];

A = [ F1; F21, F22, F23; F3 ];

f = [feval(cds.curve_func, x) + b*p; J(:,1:ndim-1)'*p; p'*J(:,ndim); p'*p-1];

% ---------------------------------------------------------

function WorkspaceInit(x,v)

% ------------------------------------------------------

function WorkspaceDone

% -------------------------------------------------------

function [res,x,v] = Adapt(x,v)
res = []; % no re-evaluations needed

%SD:continues equilibrium of odefile
