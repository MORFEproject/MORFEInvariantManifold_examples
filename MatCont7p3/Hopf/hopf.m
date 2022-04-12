function out = Hopf
%
% Hopf curve definition file for a problem in odefile
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
global hds
  [x,p,k] = rearr(arg); p = num2cell(p);
  jac = cjac(hds.func,hds.Jacobian,x,p,hds.ActiveParams);
  RED=jac*jac+k*eye(hds.nphase);
  Bord=[RED hds.borders.w;hds.borders.v' zeros(2)];
  bunit=[zeros(hds.nphase,2);eye(2)];
  vext=Bord\bunit;
  vext1=vext(hds.nphase+hds.index1(1),hds.index1(2));
  vext2=vext(hds.nphase+hds.index2(1),hds.index2(2));
  func = [feval(hds.func, 0, x, p{:}); vext1;vext2];

%--------------------------------------------------------  
function jac = jacobian(varargin)
global hds
    nap = length(hds.ActiveParams);
    xo = varargin{1}; [x,p,k] = rearr(xo); p = num2cell(p);
    J=cjac(hds.func,hds.Jacobian,x,p,hds.ActiveParams);
    RED=J*J+k*eye(hds.nphase);
    Bord=[RED hds.borders.w;hds.borders.v' zeros(2)];
    bunit=[zeros(hds.nphase,2);eye(2)];
    vext=Bord\bunit;
    wext=Bord'\bunit;
    jac=[J  cjacp(hds.func,hds.JacobianP,x,p,hds.ActiveParams) zeros(hds.nphase,1)];
    hess=chess(hds.func,hds.Jacobian,hds.Hessians,x,p,hds.ActiveParams);
    for i=1:hds.nphase
        jac(hds.nphase+1,i)=-wext(1:hds.nphase,hds.index1(1))'*(J*hess(:,:,i)+hess(:,:,i)*J)*vext(1:hds.nphase,hds.index1(2));
        jac(hds.nphase+2,i)=-wext(1:hds.nphase,hds.index2(1))'*(J*hess(:,:,i)+hess(:,:,i)*J)*vext(1:hds.nphase,hds.index2(2));
    end
    jac(hds.nphase+1,hds.nphase+nap+1)=-wext(1:hds.nphase,hds.index1(1))'*vext(1:hds.nphase,hds.index1(2));
    jac(hds.nphase+2,hds.nphase+nap+1)=-wext(1:hds.nphase,hds.index2(1))'*vext(1:hds.nphase,hds.index2(2));
    hessp=chessp(hds.func,hds.Jacobian,hds.HessiansP,x,p,hds.ActiveParams);
    for i=1:nap
        jac(hds.nphase+1,hds.nphase+i)=-wext(1:hds.nphase,hds.index1(1))'*(J*hessp(:,:,i)+hessp(:,:,i)*J)*vext(1:hds.nphase,hds.index1(2));
        jac(hds.nphase+2,hds.nphase+i)=-wext(1:hds.nphase,hds.index2(1))'*(J*hessp(:,:,i)+hessp(:,:,i)*J)*vext(1:hds.nphase,hds.index2(2));
    end    
%------------------------------------------------------
function hessians(varargin)  
%------------------------------------------------------
function varargout = defaultprocessor(varargin)
global cds hds
%elseif strcmp(arg,'defaultprocessor')
  if nargin > 2
    s = varargin{3};
    varargout{3} = s;
  end
 % compute eigenvalues?
  if (cds.options.Eigenvalues==1)
      xo = varargin{1}; [x,p] = rearr(xo); p = num2cell(p);
      d = eig(cjac(hds.func,hds.Jacobian,x,p,hds.ActiveParams));
      [Y,I] = sort(real(d));
      varargout{2} = d(I);
  else
      varargout{2}= nan;
  end  

  % all done succesfully
  varargout{1} = 0;
%-------------------------------------------------------  
function option = options
global hds
  option = contset;
  switch hds.nphase
      case 1
          option=contset(option,'IgnoreSingularity',[1 2 3 4]);
      case 2
          option=contset(option,'IgnoreSingularity',[2 3]);
      case 3
          option=contset(option,'IgnoreSingularity',3);
  end
  % Check for symbolic derivatives in odefile
  
  symord = 0; 
  if ~isempty(hds.Jacobian), symord = 1; end
  if ~isempty(hds.Hessians), symord = 2; end
  if ~isempty(hds.Der3), symord = 3; end
  if ~isempty(hds.Der4), symord = 4; end
  if ~isempty(hds.Der5), symord = 5; end

  option = contset(option, 'SymDerivative', symord);
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', [0 0 0]);

  symordp = 0;
  if ~isempty(hds.JacobianP)
      symordp = 1; 
  end
  if ~isempty(hds.HessiansP)
      symordp = 2;
  end
  option=contset(option,'SymDerivativeP',symordp);
 

% ---------------------------------------------------------------

function [out, failed] = testf(id, x, v)
global hds
[x0,p,k] = rearr(x);
p1 = num2cell(p);
jac=cjac(hds.func,hds.Jacobian,x0,p1,hds.ActiveParams);

% if any(ismember(id,[3 4]))
if any(ismember(id,[3 4]))
    RED=jac*jac+k*eye(hds.nphase);
    Bord=[RED hds.borders.w;hds.borders.v' zeros(2)];
    bunit=[zeros(hds.nphase,2);eye(2)];
    vext=Bord\bunit;
    wext=Bord'\bunit;
end

failed = [];

for i=id
  lastwarn('');
  
  switch i
  case 1 % BT      
      out(1) =k;       
  case 2 % ZH
      out(2) = det(jac);    
  case 3 % HH
      A1 = sparse(hds.BiAlt_M1_I,hds.BiAlt_M1_J,jac(hds.BiAlt_M1_V));
      A2 = sparse(hds.BiAlt_M2_I,hds.BiAlt_M2_J,jac(hds.BiAlt_M2_V));
      A3 = sparse(hds.BiAlt_M3_I,hds.BiAlt_M3_J,jac(hds.BiAlt_M3_V));
      A = A1-A2+A3;
      bigmat = [A hds.bigW; hds.bigV' hds.bigD];
      Xg = bigmat \ [zeros(hds.nphase*(hds.nphase-1)/2+1,1); 1];
      out(3) = Xg(end);
  case 4  %GH
      if k>0
          out(4) = nf_H(hds.func,hds.Jacobian,hds.Hessians,hds.Der3,x0,p1,k,vext,wext,hds.nphase,hds.ActiveParams);
      else 
          out(4) = 500;
      end
  otherwise
    error('No such testfunction');
  end 
end 

%-------------------------------------------------------
function [out, failed] = userf(userinf, id, x, v)
global hds
dim =size(id,2);
failed = [];
out(dim) = 0;
for i=1:dim
  lastwarn('');
  [x0,p] = rearr(x); p = num2cell(p);
  if (userinf(i).state==1)
      out(i)=feval(hds.user{id(i)},0,x0,p{:});
  else
      out(i)=0;
  end
  if ~isempty(lastwarn)
    failed = [failed i];
  end
end
% ---------------------------------------------------------------

function [failed,s] = process(id, x, v, s)
global cds hds
ndim = cds.ndim;
[x0,p,k] = rearr(x); p = num2cell(p);
jac=cjac(hds.func,hds.Jacobian,x0,p,hds.ActiveParams);

% WM: Removed SL array
fprintf('label = %s, x = ', s.label); printv(x);
if any(ismember(id,[4]))
    RED=jac*jac+k*eye(hds.nphase);
    Bord=[RED hds.borders.w;hds.borders.v' zeros(2)];
    bunit=[zeros(hds.nphase,2);eye(2)];
    vext=Bord\bunit;
    wext=Bord'\bunit;
end

switch id
  case 1 % BT
    s.data.c = nf_BT(hds.func,hds.Jacobian,hds.Hessians,x0,p,hds.nphase);
    fprintf('(a,b)=(%d, %d)\n',s.data.c);
    s.msg  = sprintf('Bogdanov-Takens point');
  case 2 % ZH
    s.data.c = nf_ZH(hds.func,hds.Jacobian,hds.Hessians,hds.Der3,x0,p,hds.nphase);
    if all(real(s.data.c) == 0)
        fprintf('Zero-Neutral Saddle\n');
        s.msg  = sprintf('Zero-Hopf point: neutral saddle');
    else
        fprintf('(s,theta,E0)=(%d, %d, %d)\n',real(s.data.c));
        s.msg  = sprintf('Zero-Hopf point');
    end
  case 3 % HH
    s.data.c = nf_HH(hds.func,hds.Jacobian,hds.Hessians,hds.Der3,hds.Der4,hds.Der5,x0,p,hds.nphase);
    if ~isempty(s.data.c)
        fprintf('(p11*p22,theta,delta)=(%d, %d, %d)\n (THETA,DELTA) = (%d, %d)\n',s.data.c);
        s.msg  = sprintf('Double Hopf point');
    else
        s.msg = sprintf('Neutral saddle');
    end
  case 4 %GH
    s.data.c = nf_GH(hds.func,hds.Jacobian,hds.Hessians,hds.Der3,hds.Der4,hds.Der5,x0,p,k,vext,wext,hds.nphase);
    fprintf('l2=%d\n',s.data.c);
    s.msg  = sprintf('Generalized Hopf point');  
end

% Compute eigenvalues for every singularity
J=cjac(hds.func,hds.Jacobian,x0,p,hds.ActiveParams);
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
%-------------------------------------------------------
function [S,L] = singmat
 
% 0: testfunction must vanish
% 1: testfunction must not vanish
% everything else: ignore this testfunction

  S = [  0 0 8 8
         1 0 8 8
         8 8 0 8
         8 8 8 0];
     
  L = [ 'BT'; 'ZH'; 'HH';'GH'];

%-------------------------------------------------------
function locate(id, varargin)
error('No locator defined for singularity %d', id);
  
%-------------------------------------------------------
function varargout = init(varargin)
  x = varargin{1};
  v = varargin{2};
  WorkspaceInit(x,v);

  % all done succesfully
  varargout{1} = 0;
%-------------------------------------------------------
function done()
  WorkspaceDone;
  
%------------------------------------------------------
function [res,x,v] = adapt(x,v)
global hds
nap=length(hds.ActiveParams);
[x0,p,k] = rearr(x); p1 = num2cell(p);
jac=cjac(hds.func,hds.Jacobian,x0,p1,hds.ActiveParams);
RED=jac*jac+k*eye(hds.nphase);
Bord=[RED hds.borders.w;hds.borders.v' zeros(2)];
bunit=[zeros(hds.nphase,2);eye(2)];
vext=Bord\bunit;
[vext,r]=qr(vext);
vext=vext(:,1:2);hds.borders.v=vext(1:hds.nphase,:);
wext=Bord'\bunit;
[wext,r]=qr(wext);
wext=wext(:,1:2);hds.borders.w=wext(1:hds.nphase,:);
A=[jac cjacp(hds.func,hds.JacobianP,x0,p1,hds.ActiveParams) zeros(hds.nphase,1)];
[Q,R]=qr(A');
hess=chess(hds.func,hds.Jacobian,hds.Hessians,x0,p1,hds.ActiveParams);
gx(4,hds.nphase) = 0;
for i=1:hds.nphase
    gx(1,i)=-wext(1:hds.nphase,1)'*(jac*hess(:,:,i)+hess(:,:,i)*jac)*vext(1:hds.nphase,1);
    gx(2,i)=-wext(1:hds.nphase,1)'*(jac*hess(:,:,i)+hess(:,:,i)*jac)*vext(1:hds.nphase,2);
    gx(3,i)=-wext(1:hds.nphase,2)'*(jac*hess(:,:,i)+hess(:,:,i)*jac)*vext(1:hds.nphase,1);
    gx(4,i)=-wext(1:hds.nphase,2)'*(jac*hess(:,:,i)+hess(:,:,i)*jac)*vext(1:hds.nphase,2);
end
gk(1,1)=-wext(1:hds.nphase,1)'*vext(1:hds.nphase,1);
gk(2,1)=-wext(1:hds.nphase,1)'*vext(1:hds.nphase,2);
gk(3,1)=-wext(1:hds.nphase,2)'*vext(1:hds.nphase,1);
gk(4,1)=-wext(1:hds.nphase,2)'*vext(1:hds.nphase,2);
hessp=chessp(hds.func,hds.Jacobian,hds.HessiansP,x0,p1,hds.ActiveParams);
gp(4,nap) = 0;
for i=1:nap
    gp(1,i)=-wext(1:hds.nphase,1)'*(jac*hessp(:,:,i)+hessp(:,:,i)*jac)*vext(1:hds.nphase,1);
    gp(2,i)=-wext(1:hds.nphase,1)'*(jac*hessp(:,:,i)+hessp(:,:,i)*jac)*vext(1:hds.nphase,2);
    gp(3,i)=-wext(1:hds.nphase,2)'*(jac*hessp(:,:,i)+hessp(:,:,i)*jac)*vext(1:hds.nphase,1);
    gp(4,i)=-wext(1:hds.nphase,2)'*(jac*hessp(:,:,i)+hessp(:,:,i)*jac)*vext(1:hds.nphase,2);
end
A=[A;gx gp gk]*Q;
Jres=A(1+hds.nphase:end,1+hds.nphase:end)';
[Q,R,E]=qr(Jres);
E=E(:,1:2)';
index=[1 1;1 2;2 1;2 2];
[I,J]=find(E);
hds.index1=index(J(1),:);
hds.index2=index(J(2),:);

if hds.nphase > 3
    A1 = sparse(hds.BiAlt_M1_I,hds.BiAlt_M1_J,jac(hds.BiAlt_M1_V));
    A2 = sparse(hds.BiAlt_M2_I,hds.BiAlt_M2_J,jac(hds.BiAlt_M2_V));
    A3 = sparse(hds.BiAlt_M3_I,hds.BiAlt_M3_J,jac(hds.BiAlt_M3_V));
    A = A1-A2+A3;
    bigmat = [A hds.bigW; hds.bigV' hds.bigD];
    if hds.adaptchoice == 2
        bigmat = bigmat';
    end
    XY = bigmat \ [zeros(hds.nphase*(hds.nphase-1)/2,2); eye(2,2)];
    yopt1 = [XY(end,2); -XY(end,1)];
    yopt2 = [XY(end-1,2); -XY(end-1,1)];
    if norm(yopt1) > norm(yopt2)
        smally = yopt1;
    else
        smally = yopt2;
    end
    X1tmp = XY*smally;
    X1tmp = X1tmp / norm(X1tmp);
    X1 = X1tmp(1:end-2,1);
    alpha = -X1tmp' * XY(:,2);
    beta = X1tmp' * XY(:,1);
    X2tmp = alpha * XY(:,1) + beta * XY(:,2);
    X2tmp = X2tmp / norm(X2tmp);
    hds.bigD = [zeros(1,2); X2tmp(end-1:end)'];
    if hds.adaptchoice == 1
        hds.bigV = [X1 X2tmp(1:end-2)];
        hds.adaptchoice = 2;
    else
        hds.bigW = [X1 X2tmp(1:end-2)];
        hds.adaptchoice = 1;
        hds.bigD = hds.bigD';
    end
end

res = 1; % no re-evaluations needed

 

%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------

function [x,p,k] = rearr(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x) and parameters (p)
global hds
p = hds.P0;
p(hds.ActiveParams) = x0((hds.nphase+1):(end-1));
x = x0(1:hds.nphase);
k=x0(end);
   
% ---------------------------------------------------------

function WorkspaceInit(x,v)

% ------------------------------------------------------

function WorkspaceDone

% -------------------------------------------------------

%SD:continues Hopf of odefile
