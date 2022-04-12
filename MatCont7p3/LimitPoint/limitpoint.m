function out = limitpoint
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

%----------------------------------------------------
function func = curve_func(arg)
global lpds
  [x,p] = rearr(arg); p = num2cell(p);
  Bord=[cjac(lpds.func,lpds.Jacobian,x,p,lpds.ActiveParams) lpds.borders.w;lpds.borders.v' 0];
  bunit=[zeros(lpds.nphase,1);1];
  vext=Bord\bunit;
  func = [feval(lpds.func, 0, x, p{:}) ; vext(end)];
  
%---------------------------------------------------     
function jac = jacobian(varargin)
global lpds
  nap = length(lpds.ActiveParams);
  xo = varargin{1}; [x,p] = rearr(xo); p = num2cell(p);
  jac=cjac(lpds.func,lpds.Jacobian,x,p,lpds.ActiveParams);
  Bord=[jac lpds.borders.w;lpds.borders.v' 0];
  bunit=[zeros(lpds.nphase,1);1];
  vext=Bord\bunit;
  wext=Bord'\bunit;
  jac = [jac cjacp(lpds.func,lpds.JacobianP,x,p,lpds.ActiveParams)]; 
  hess=chess(lpds.func,lpds.Jacobian,lpds.Hessians,x,p,lpds.ActiveParams);
  for i=1:lpds.nphase
      jac(lpds.nphase+1,i)=-wext(1:lpds.nphase)'*hess(:,:,i)*vext(1:lpds.nphase);
  end
  hessp=chessp(lpds.func,lpds.Jacobian,lpds.HessiansP,x,p,lpds.ActiveParams);
  for i=1:nap
      jac(lpds.nphase+1,lpds.nphase+i)=-wext(1:lpds.nphase)'*hessp(:,:,i)*vext(1:lpds.nphase);
  end
%---------------------------------------------------
function hessians(varargin)
%---------------------------------------------------
function varargout = defaultprocessor(varargin)
global  cds lpds
  if nargin > 2
    s = varargin{3};
    varargout{3} = s;
  end
 % compute eigenvalues?
  if (cds.options.Eigenvalues==1)
      xo = varargin{1}; [x,p] = rearr(xo); p = num2cell(p);
      d = eig(cjac(lpds.func,lpds.Jacobian,x,p,lpds.ActiveParams));
      [Y,I] = sort(real(d));
      varargout{2} = d(I);
  else
      varargout{2} = nan;
  end  

  % all done succesfully
  varargout{1} = 0;
%----------------------------------------------------  
function option = options
global lpds 

  bp = size(lpds.BranchParams,2);
  option = contset;
  switch lpds.nphase
      case 1
          option=contset(option,'IgnoreSingularity',[bp+1 bp+2]);
      case 2
          option=contset(option,'IgnoreSingularity',[bp+2]);
  end

  % Check for symbolic derivatives in odefile
  
  symord = 0; 
  if ~isempty(lpds.Jacobian), symord = 1; end
  if ~isempty(lpds.Hessians), symord = 2; end
  if ~isempty(lpds.Der3), symord = 3; end
  if ~isempty(lpds.Der4), symord = 4; end
  if ~isempty(lpds.Der5), symord = 5; end

  option = contset(option, 'SymDerivative', symord);
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', [0 0 0]);

  symordp = 0;
  if ~isempty(lpds.JacobianP)
      symordp = 1;
  end
  if ~isempty(lpds.HessiansP)
      symordp = 2;
  end
  option = contset(option,'SymDerivativeP',symordp);

% -------------------------------------------------------
function [out, failed] = testf(id, x, v)
global lpds
  [x0,p] = rearr(x); p = num2cell(p);
  jac=cjac(lpds.func,lpds.Jacobian,x0,p,lpds.ActiveParams);
  Bord=[jac lpds.borders.w;lpds.borders.v' 0];
  bunit=[zeros(lpds.nphase,1);1];
  vext=Bord\bunit;
  vext=vext(1:lpds.nphase);
  wext=Bord'\bunit;
  wext=wext(1:lpds.nphase);
  failed = [];

if lpds.nphase > 1
if ismember(1,id)% case 1 % BT
    out(1) = wext'*vext;
end
end
if lpds.nphase > 2
if ismember(2,id)% case 2 % ZH
    A = jac;
    A1=sparse(lpds.BiAlt_M1_I,lpds.BiAlt_M1_J,A(lpds.BiAlt_M1_V));
    A2=sparse(lpds.BiAlt_M2_I,lpds.BiAlt_M2_J,A(lpds.BiAlt_M2_V));
    A3=sparse(lpds.BiAlt_M3_I,lpds.BiAlt_M3_J,A(lpds.BiAlt_M3_V));
    A = A1-A2+A3;
    bigmat = [A lpds.bigW; lpds.bigV' lpds.bigD];
    
    Xg = bigmat \ [zeros(lpds.nphase*(lpds.nphase-1)/2,1); 1];
    out(2) = Xg(end);
end
end
if ismember(3,id)% case 3 % CP
     out(3) = nf_LP(lpds.func,lpds.Jacobian,lpds.Hessians,x0,p,vext/norm(vext),wext/(wext'*vext),lpds.nphase);
end
if id(id>3)%  %BP
    bp = size(lpds.BranchParams,2);
    out(3+(1:bp))=wext'*cjacbr(lpds.func,lpds.JacobianP,x0,p,lpds.ActiveParams,lpds.BranchParams);
end
%------------------------------------------------------
function [out, failed] = userf(userinf, id, x, v)
global lpds
dim =size(id,2);
failed = [];
out(dim) = 0;
for i=1:dim
  lastwarn('');
  [x0,p] = rearr(x); p = num2cell(p);
  if (userinf(i).state==1)
      out(i)=feval(lpds.user{id(i)},0,x0,p{:});
  else
      out(i)=0;
  end
  if ~isempty(lastwarn)
    failed = [failed i];
  end
end
%---------------------------------------------------------
function [failed,s] = process(id, x, v, s)
global cds lpds
  ndim = cds.ndim;
  [x0,p] = rearr(x); p = num2cell(p);
% WM: Removed SL array
fprintf('label = %s, x = ', s.label); printv(x);

switch id
  case 1 % BT
    s.data.c = nf_BT(lpds.func,lpds.Jacobian,lpds.Hessians,x0,p,lpds.nphase);
    fprintf('(a,b)=(%d, %d)\n',s.data.c);
    s.msg  = sprintf('Bogdanov-Takens point');
  case 2 % ZH
    s.data.c = nf_ZH(lpds.func,lpds.Jacobian,lpds.Hessians,lpds.Der3,x0,p,lpds.nphase);
    if ~all(real(s.data.c) == 0)
        fprintf('(s,theta,E0)=(%d, %d, %d)\n',real(s.data.c));
        s.msg = sprintf('Zero-Hopf point');
    else
        fprintf('Zero-Neutral Saddle\n');
        s.msg  = sprintf('Zero-Hopf point: neutral saddle');
    end
  case 3 % CP
    s.data.c = nf_CP(lpds.func,lpds.Jacobian,lpds.Hessians,lpds.Der3,x0,p,lpds.nphase);
    fprintf('c=%d\n',s.data.c);
    s.msg  = sprintf('Cusp point');
  otherwise
    s.msg = sprintf('Branch point (parameter %d)',lpds.BranchParams(id-3));
end

% Compute eigenvalues for every singularity
[x0,p] = rearr(x); p = num2cell(p); 
J=cjac(lpds.func,lpds.Jacobian,x0,p,lpds.ActiveParams); 
if ~issparse(J)
  [v,d]=eig(J);
else
  opt.disp=0;
  % WM: fixed a bug (incompatability between MatLab 6.0 and 5.5?)
  [v,d]=eigs(J,min(6,ndim-2),'lm',opt);
end

s.data.evec = v;
s.data.eval = diag(d)';

failed = 0;
%--------------------------------------------------------
function  [S,L] = singmat    
global lpds 
% 0: testfunction must vanish
% 1: testfunction must not vanish
% everything else: ignore this testfunction

  S = [  0 0 8 
         1 0 8 
         8 8 0 ];
     if (lpds.nphase==2)
         S(1,2)=8;
     end
  bp = size(lpds.BranchParams,2);                    
  for i = 1:bp
      S(3+i,end+1)=0;
      S(3+i,1:end-1)=8;
      S(1:end-1,3+i)=8;
  end
%  L = [ 'BT '; 'ZH '; 'CP '];
%  for i = 1:bp
%      L(3+i,:)= strcat('BP',num2str(lpds.BranchParams(i)));
%  end
  L = { 'BT'; 'ZH'; 'CP'};
  for i = 1:bp
      L{3+i}= strcat('BP',num2str(lpds.BranchParams(i)));
  end
  % recast the cell array back to a character array
  L = char ( L );
%------------------------------------------------------
function [x,v] = locate(id, varargin)
error('No locator defined for singularity %d', id);
%------------------------------------------------------
function varargout = init(varargin)
  x = varargin{1};
  v = varargin{2};
  WorkspaceInit(x,v);

  % all done succesfully
  varargout{1} = 0;
%--------------------------------------------------------
function done
  WorkspaceDone;
%---------------------------------------------------------
function [res,x,v] = adapt(x,v)
global lpds 
[x1,p] =rearr(x); p = num2cell(p);
jac = cjac(lpds.func,lpds.Jacobian,x1,p,lpds.ActiveParams);
Bord=[jac lpds.borders.w;lpds.borders.v' 0];
bunit=[zeros(lpds.nphase,1);1];
vext=Bord\bunit;
wext=Bord'\bunit;
%ERROR OR WARNING
lpds.borders.v=vext(1:lpds.nphase)/norm(vext(1:lpds.nphase));
lpds.borders.w=wext(1:lpds.nphase)/norm(wext(1:lpds.nphase));

if lpds.nphase > 2
    A1 = sparse(lpds.BiAlt_M1_I,lpds.BiAlt_M1_J,jac(lpds.BiAlt_M1_V));
    A2 = sparse(lpds.BiAlt_M2_I,lpds.BiAlt_M2_J,jac(lpds.BiAlt_M2_V));
    A3 = sparse(lpds.BiAlt_M3_I,lpds.BiAlt_M3_J,jac(lpds.BiAlt_M3_V));
    A = A1-A2+A3;
    [Q,R,E] = qr(full(A));
    if lpds.bigW' * Q(:,end) < 0
        lpds.bigW = -Q(:,end);
    else
        lpds.bigW = Q(:,end);
    end
    if lpds.bigV' * E(:,end) < 0
        lpds.bigV = -E(:,end);
    else
        lpds.bigV = E(:,end);
    end
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
global lpds
p = lpds.P0;
p(lpds.ActiveParams) = x0((lpds.nphase+1):end);
x = x0(1:lpds.nphase);

    
% ---------------------------------------------------------

function WorkspaceInit(x,v)
global lpds 

% calculate some matrices to efficiently compute bialternate products
% (without loops)
n = lpds.nphase;
a = reshape(1:(n^2),n,n);
[bia,bin,bip] = bialt(a);
if any(any(bip))
    [lpds.BiAlt_M1_I,lpds.BiAlt_M1_J,lpds.BiAlt_M1_V] = find(bip);
else
    lpds.BiAlt_M1_I=1;lpds.BiAlt_M1_J=1;lpds.BiAlt_M1_V=n^2+1;
end    
if any(any(bin))
    [lpds.BiAlt_M2_I,lpds.BiAlt_M2_J,lpds.BiAlt_M2_V] = find(bin);
else
     lpds.BiAlt_M2_I=1;lpds.BiAlt_M2_J=1;lpds.BiAlt_M2_V=n^2+1;
end
if any(any(bia))
    [lpds.BiAlt_M3_I,lpds.BiAlt_M3_J,lpds.BiAlt_M3_V] = find(bia);
else
    lpds.BiAlt_M3_I=1;lpds.BiAlt_M3_J=1;lpds.BiAlt_M3_V=n^2+1;
end

% ------------------------------------------------------

function WorkspaceDone

% -------------------------------------------------------
%SD:continues equilibrium of odefile
