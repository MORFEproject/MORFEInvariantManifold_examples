function out = equilibrium
%
% Equilibrium curve definition file for a problem in odefile
% 

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
%-------------------------------------------------------
function func = curve_func(arg)
global eds 
  [x,p] = rearr(arg);  p = num2cell(p);
  func = feval(eds.func,0, x, p{:});  
%---------------------------------------------------------------
function jac = jacobian(varargin)
global eds
    [x,p] = rearr(varargin{1});    p = num2cell(p);
    jac = [cjac(eds.func,eds.Jacobian,x,p,eds.ActiveParams) cjacp(eds.func,eds.JacobianP,x,p,eds.ActiveParams)];
%---------------------------------------------------------------    
function hess = hessians(varargin)  
global cds eds
    xo = varargin{1}; [x,p] =  rearr(xo);p=num2cell(p);
    hh = chess(eds.func,eds.Jacobian,eds.Hessians,x,p,eds.ActiveParams);
    hp = chessp(eds.func,eds.Jacobian,eds.HessiansP,x,p,eds.ActiveParams);
    x1 = xo; x1(cds.ndim) = x1(cds.ndim) - cds.options.Increment;
    x2 = xo; x2(cds.ndim) = x2(cds.ndim) + cds.options.Increment;
    hpp = (cjac(cds.curve_func,cds.curve_jacobian,x2,[]) - cjac(cds.curve_func,cds.curve_jacobian,x1,[])) / (2*cds.options.Increment);
    for i = 1:cds.ndim-1
        hess(:,:,i) = [ hh(:,:,i) hpp(:,i)];
    end
    hess(:,:,cds.ndim) = [ hp(:,:) hpp(:,cds.ndim)]; 
%---------------------------------------------------------------
function varargout = defaultprocessor(varargin)
global eds cds
  if nargin > 2
    s = varargin{3};
    s.data.v = eds.v;
    varargout{3} = s;
  end

  % compute eigenvalues?
  if (cds.options.Eigenvalues==1)
      xo = varargin{1}; [x,p] = rearr(xo); p = num2cell(p);
      d = eig(cjac(eds.func,eds.Jacobian,x,p,eds.ActiveParams));
      [Y,I] = sort(real(d));
      varargout{2} = d(I);
  else
      varargout{2}=nan;
  end  
  % all done succesfully
  varargout{1} = 0;
%-------------------------------------------------------------
function option = options
global eds cds
  option = contset;
  if cds.ndim<3
          option=contset(option,'IgnoreSingularity',2);
  end
  % Check for symbolic derivatives in odefile
  
  symord = 0; 
  if ~isempty(eds.Jacobian), symord = 1; end
  if ~isempty(eds.Hessians), symord = 2; end
  if ~isempty(eds.Der3), symord = 3; end
  if ~isempty(eds.Der4), symord = 4; end
  if ~isempty(eds.Der5), symord = 5; end

  option = contset(option, 'SymDerivative', symord);
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', [1 0 0]);

  symordp = 0;
  if ~isempty(eds.JacobianP), symordp = 1; end
  if ~isempty(eds.HessiansP),  symordp = 2; end
  option = contset(option, 'SymDerivativeP', symordp);

%----------------------------------------------------------------
function [out, failed] = testf(id, x, v)
global cds eds
ndim = cds.ndim;

if any(ismember(id,[1 2]))
    J = cjac(cds.curve_func,cds.curve_jacobian,x,[]);
end


out(3) = 0;
failed = [];

for i=id
  lastwarn('');
  
  switch i
  case 1 % BP
    % Jacobian extended with bordering vectors v and w
    B = [J; v'];
    out(1) = det(B);
    
  case 2 % H
    A=J(:,1:ndim-1);
    A(:,end+1)=0;
    A1=sparse(eds.BiAlt_M1_I,eds.BiAlt_M1_J,A(eds.BiAlt_M1_V));
    A2=sparse(eds.BiAlt_M2_I,eds.BiAlt_M2_J,A(eds.BiAlt_M2_V));
    A3=sparse(eds.BiAlt_M3_I,eds.BiAlt_M3_J,A(eds.BiAlt_M3_V));
    A = A1-A2+A3;
    bigmat = [A eds.bigW; eds.bigV' eds.bigD];
    
    Xg = bigmat \ [zeros(eds.nphase*(eds.nphase-1)/2,1); 1];
    out(2) = Xg(end);
    
  case 3 % LP
    out(3) = v(end);
    
  otherwise
    error('No such testfunction');
  end
  if ~isempty(lastwarn)
    failed = [failed id];
  end
  
end
%-----------------------------------------------------------------
function [out, failed] = userf(userinf, id, x, v)
global eds
dim =size(id,2);
failed = [];
out(dim) = 0;
for i=1:dim
  lastwarn('');
  [x0,p] = rearr(x); p = num2cell(p);
  if (userinf(i).state==1)
      out(i)=feval(eds.user{id(i)},0,x0,p{:});
  else
      out(i)=0;
  end
  if ~isempty(lastwarn)
    failed = [failed i];
  end
end
%---------------------------------------------------------------------
function [failed,s] = process(id, x, v, s)
global cds
ndim = cds.ndim;

% WM: Removed SL array
% fprintf('label = %s, x = ', s.label); printv(x);

switch id
  case 1 % BP
    s.data.v = v;
    s.msg  = sprintf('Branch point');
    fprintf('label = %s, x = ', s.label); printv(x);
  case 2 % H
    s.data.lyapunov = lyapunov(x);
    if strcmp(s.data.lyapunov,'Neutral saddle')
        s.msg  = sprintf('Neutral Saddle Equilibrium');
        fprintf('label = %s, x = ', s.label); printv(x);
    else
        s.msg  = sprintf('Hopf');
        fprintf('label = %s, x = ', s.label); printv(x);
        fprintf('First Lyapunov coefficient = %d\n', s.data.lyapunov);
    end
  case 3 % LP
    s.data.a = a_lp(x);
    s.msg  = sprintf('Limit point');
    fprintf('label = %s, x = ', s.label); printv(x);
    fprintf('a=%d\n',s.data.a);
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

  S = [  0 8 8
         8 0 8
         1 8 0 ];

  L = [ 'BP'; 'H '; 'LP' ];


  %elseif strcmp(arg, 'locate')
%--------------------------------------------------------
function [x,v] = locate(id, x1, v1, x2, v2)
switch id
  case 1
    [x,v] = locateBP(id, x1, v1, x2, v2);
  otherwise
    error('No locator defined for singularity %d', id);
end
%---------------------------------------------------------
function varargout = init(varargin)
  x = varargin{1};
  v = varargin{2};
  WorkspaceInit(x,v);

  % all done succesfully
  varargout{1} = 0;
%---------------------------------------------------------
function done
  WorkspaceDone;

%----------------------------------------------------------  
function [res,x,v] = adapt(x,v)
global eds cds
jac = cjac(cds.curve_func,cds.curve_jacobian,x,[]);
jac=jac(:,1:cds.ndim-1);
jac(:,end+1)=0;
A1 = sparse(eds.BiAlt_M1_I,eds.BiAlt_M1_J,jac(eds.BiAlt_M1_V));
A2 = sparse(eds.BiAlt_M2_I,eds.BiAlt_M2_J,jac(eds.BiAlt_M2_V));
A3 = sparse(eds.BiAlt_M3_I,eds.BiAlt_M3_J,jac(eds.BiAlt_M3_V));
A = A1-A2+A3;
[Q,R,E] = qr(full(A));
if eds.bigW' * Q(:,end) < 0
    eds.bigW = -Q(:,end);
else
    eds.bigW = Q(:,end);
end
if eds.bigV' * E(:,end) < 0
    eds.bigV = -E(:,end);
else
    eds.bigV = E(:,end);
end

res = 1; 

  


%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------

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
%if ~issparse(J)
%  [v,d]=eig(J(:,1:ndim-1)');
%else
%  opt.disp=0;
%  [v,d]=eigs(J(:,1:ndim-1)', 'SM', opt);
%end
%[y,i]=min(abs(diag(d)));
%p = v(:,i);
[U,S,V]=svd(full(J)',0);
p=V(:,end);
b = 0;
x = 0.5*(x1+x2);
v = 0.5*(v1+v2);
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
  if norm(du) < 100*cds.options.VarTolerance & norm(f) < cds.options.FunTolerance 
      return; 
  end
  i = i+1;
end
x=[];


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
global cds eds

% calculate some matrices to efficiently compute bialternate products (without loops)
n = cds.ndim-1;
a = reshape(1:(n^2),n,n);
[bia,bin,bip] = bialt(a);
if any(any(bip))
    [eds.BiAlt_M1_I,eds.BiAlt_M1_J,eds.BiAlt_M1_V] = find(bip);
else
    eds.BiAlt_M1_I=1;eds.BiAlt_M1_J=1;eds.BiAlt_M1_V=n^2+1;
end    
if any(any(bin))
    [eds.BiAlt_M2_I,eds.BiAlt_M2_J,eds.BiAlt_M2_V] = find(bin);
else
     eds.BiAlt_M2_I=1;eds.BiAlt_M2_J=1;eds.BiAlt_M2_V=n^2+1;
end
if any(any(bia))
    [eds.BiAlt_M3_I,eds.BiAlt_M3_J,eds.BiAlt_M3_V] = find(bia);
else
    eds.BiAlt_M3_I=1;eds.BiAlt_M3_J=1;eds.BiAlt_M3_V=n^2+1;
end

% ------------------------------------------------------

function WorkspaceDone

% -------------------------------------------------------


%SD:continues equilibrium of odefile
