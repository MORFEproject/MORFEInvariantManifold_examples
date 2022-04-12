function out = heteroclinic
%
% 
global hetds cds
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
global hetds

[x,x0,x1,p,T,eps0,eps1,YS,YU] = rearr(arg);
func = BVP_Het(x,x0,x1,p,T,eps0,eps1,YS,YU);  

%------------------------------------------------------
function varargout = jacobian(varargin)
global hetds BigJac

[x,x0,x1,p,T,eps0,eps1,YS,YU] = rearr(varargin{1});
varargout{1} = BVP_Het_jac(hetds.func,x,x0,x1,p,T,eps0,eps1,YS,YU);
BigJac = varargout{1};

%-----------------------------------------------------

function hessians(varargin)

%------------------------------------------------------

function varargout = defaultprocessor(varargin)
global hetds cds

[x,x0,x1,p,T,eps0,eps1,YS,YU] = rearr(varargin{1});
v = rearr(varargin{2});

hetds.ndim = length(varargin{1});

% update
if ~isempty(hetds.ups)
    hetds.upold = hetds.ups;
end
hetds.ups = reshape(x,hetds.nphase,hetds.tps);
hetds.vps = reshape(v,hetds.nphase,hetds.tps);

p1 = num2cell(p);
for i=1:hetds.tps      
    hetds.upoldp(:,i) = 2*T*feval(hetds.func, 0, hetds.ups(:,i),p1{:});
end
hetds.T = T;
hetds.eps0 = eps0;
hetds.eps1 = eps1;
hetds.YS = YS;
hetds.YU = YU;
hetds.x0 = x0;
hetds.x1 = x1;

% Update dimensions
% -----------------
A1 = cjac(hetds.func,hetds.Jacobian,hetds.x0,p1,hetds.ActiveParams);
D1 = eig(A1);
A2 = cjac(hetds.func,hetds.Jacobian,hetds.x1,p1,hetds.ActiveParams);
D2 = eig(A2);

if nargin > 2
    % set data in special point structure
    s = varargin{3};
    s.data.timemesh = hetds.msh;
    s.data.ntst = hetds.ntst;
    s.data.ncol = hetds.ncol;
    s.data.parametervalues = p1;
    s.data.T = T;
    varargout{3} = s;
end

% all done succesfully
varargout{1} = 0;
varargout{2} = hetds.msh';
  
if (cds.options.Eigenvalues==1)
    varargout{2} = [varargout{2}; D1; D2];
end



%-------------------------------------------------------
function option = options
global hetds 

% Check for symbolic derivatives in odefile
symord = 0; 
if ~isempty(hetds.Jacobian), symord = 1; end
if ~isempty(hetds.Hessians), symord = 2; end
if ~isempty(hetds.Der3), symord = 3; end
if  ~isempty(hetds.Der4), symord = 4; end
if ~isempty(hetds.Der5), symord = 5; end

option = contset;
option = contset(option, 'SymDerivative', symord);
option = contset(option, 'Workspace', 1);
option = contset(option, 'Locators', zeros(1,13));
symordp = 0;
if  ~isempty(hetds.JacobianP), symordp = 1; end
if ~isempty(hetds.HessiansP),  symordp = 2; end
option = contset(option, 'SymDerivativeP', symordp);

%------------------------------------------------------  
function [out, failed] = testf(id, x0, v)

out = [];
failed = [];


%-------------------------------------------------------------
function [out, failed] = userf( userinf, id, x, v)

global hetds
dim =size(id,2);
failed = [];
[x,x0,x1,p,T,eps0,eps1,YS,YU] = rearr(x); p = num2cell(p);
out(dim) = 0;
for i=1:dim
  lastwarn('');
  if (userinf(i).state==1)
      out(i)=feval(hetds.user{id(i)},0,x0,p{:});
  else
      out(i)=0;
  end
  if ~isempty(lastwarn)
    failed = [failed i];
  end
end
%-----------------------------------------------------------------
function [failed,s] = process(id, x, v, s)

failed = 0;   
   
%-------------------------------------------------------------  
function [S,L] = singmat
global hetds

S = 0;
L = [];

%--------------------------------------------------------
function [x,v] = locate(id, varargin)
error('No locator defined for singularity %d', id);

%----------------------------------------------------------
function varargout = init(varargin)
WorkspaceInit(varargin{1:2});
% all done succesfully
varargout{1} = 0;
  
%-----------------------------------------------------------
function done

%-----------------------------------------------------------
function [res,x,v] = adapt(x,v)
global cds hetds

cds.adapted = 1;
YU = hetds.YU;
YS = hetds.YS;

[x,v] = Het_adapt_mesh(x,v);

Q0S = hetds.oldStableQ;
QbS1 = Q0S(:,1:hetds.nneg);
QbS2 = Q0S(:,hetds.nneg+1:end);
Q0U = hetds.oldUnstableQ;
QbU1 = Q0U(:,1:hetds.npos);
QbU2 = Q0U(:,hetds.npos+1:end);

[Q1S,S1,R1] = svd(QbS1 + QbS2*YS);
[Q1U,S1,R1] = svd(QbU1 + QbU2*YU);

%stable case
Q1S1=Q1S(:,1:hetds.nneg);
Q1S2=Q1S(:,hetds.nneg+1:end);

Q1S11=Q1S1*((Q1S1'*Q1S1)\(Q1S1'*QbS1));
for j=1:hetds.nneg
    Q1S11(:,j)=Q1S11(:,j)/norm(Q1S11(:,j));
    for i=j+1:hetds.nneg
        Q1S11(:,i)=Q1S11(:,i)-Q1S11(:,j)*(Q1S11(:,j)'*Q1S11(:,i));
    end
end
Q1S21=Q1S2*((Q1S2'*Q1S2)\(Q1S2'*QbS2));
for j=1:hetds.nphase-hetds.nneg
    Q1S21(:,j)=Q1S21(:,j)/norm(Q1S21(:,j));
    for i=j+1:hetds.nphase-hetds.nneg
        Q1S21(:,i)=Q1S21(:,i)-Q1S21(:,j)*(Q1S21(:,j)'*Q1S21(:,i));
    end
end
Q1S=[Q1S11 Q1S21];

%unstable case
Q1U1=Q1U(:,1:hetds.npos);
Q1U2=Q1U(:,hetds.npos+1:end);

Q1U11=Q1U1*((Q1U1'*Q1U1)\(Q1U1'*QbU1));
for j=1:hetds.npos
    Q1U11(:,j)=Q1U11(:,j)/norm(Q1U11(:,j));
    for i=j+1:hetds.npos
        Q1U11(:,i)=Q1U11(:,i)-Q1U11(:,j)*(Q1U11(:,j)'*Q1U11(:,i));
    end
end
Q1U21=Q1U2*((Q1U2'*Q1U2)\(Q1U2'*QbU2));
for j=1:hetds.nphase-hetds.npos
    Q1U21(:,j)=Q1U21(:,j)/norm(Q1U21(:,j));
    for i=j+1:hetds.nphase-hetds.npos
        Q1U21(:,i)=Q1U21(:,i)-Q1U21(:,j)*(Q1U21(:,j)'*Q1U21(:,i));
    end
end
Q1U=[Q1U11 Q1U21];


hetds.oldStableQ = Q1S;
hetds.oldUnstableQ = Q1U;

hetds.YS = zeros(size(hetds.YS));
hetds.YU = zeros(size(hetds.YU));

x(end-size(hetds.YS,1)*size(hetds.YS,2)-size(hetds.YU,1)*size(hetds.YU,2)+1:end)=zeros(size(hetds.YS,1)*size(hetds.YS,2)+size(hetds.YU,1)*size(hetds.YU,2),1);

res = 1;
% ---------------------------------------------------------------

%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------


%--------------------------------------------------------------
function [x,x0,x1,p,T,eps0,eps1,YS,YU] = rearr(y)
% Rearranges x1 into all of its components
global hetds

x = y(hetds.coords);
x0 = y(hetds.ncoords+1:hetds.ncoords+hetds.nphase);
x1 = y(hetds.ncoords+hetds.nphase+1:hetds.ncoords+2*hetds.nphase);

p = hetds.P0;
p(hetds.ActiveParams) = y(hetds.PeriodIdx+1:hetds.PeriodIdx+length(hetds.ActiveParams));
idx = hetds.PeriodIdx+length(hetds.ActiveParams)+1;

if hetds.extravec(1)
    T = y(idx);
    idx = idx+1;
else
    T = hetds.T;
end

if hetds.extravec(2)
    eps0 = y(idx);
    idx = idx+1;
else
    eps0 = hetds.eps0;
end

if hetds.extravec(3)
    eps1 = y(idx);
    idx = idx+1;
else
    eps1 = hetds.eps1;
end

YU = reshape(y(idx:idx+(hetds.nphase-hetds.npos)*hetds.npos-1),hetds.nphase-hetds.npos,hetds.npos);
idx = idx + (hetds.nphase-hetds.npos)*hetds.npos;
YS = reshape(y(idx:idx+(hetds.nphase-hetds.nneg)*hetds.nneg-1),hetds.nphase-hetds.nneg,hetds.nneg);
    

% -------------------------------------------------------------



% ---------------------------------------------------------------
function WorkspaceInit(x,v)
global hetds

hetds.cols_p1 = 1:(hetds.ncol+1);
hetds.cols_p1_coords = 1:(hetds.ncol+1)*hetds.nphase;
hetds.ncol_coord = hetds.ncol*hetds.nphase;
hetds.col_coords = 1:hetds.ncol*hetds.nphase;
hetds.coords = 1:hetds.ncoords;
hetds.pars = hetds.ncoords+2*hetds.nphase+(1:length(hetds.ActiveParams));
hetds.tsts = 1:hetds.ntst;
hetds.cols = 1:hetds.ncol;
hetds.phases = 1:hetds.nphase;
hetds.ntstcol = hetds.ntst*hetds.ncol;

hetds.idxmat = reshape(fix((1:((hetds.ncol+1)*hetds.ntst))/(1+1/hetds.ncol))+1,hetds.ncol+1,hetds.ntst);
hetds.dt = hetds.msh(hetds.tsts+1)-hetds.msh(hetds.tsts);

hetds.wp = kron(hetds.wpvec',eye(hetds.nphase));
hetds.pwwt = kron(hetds.wt',eye(hetds.nphase));
hetds.pwi = hetds.wi(ones(1,hetds.nphase),:);

hetds.wi = nc_weight(hetds.ncol)';

[hetds.bialt_M1,hetds.bialt_M2,hetds.bialt_M3,hetds.bialt_M4]=bialtaa(hetds.nphase);


% ------------------------------------------------------

function [x,v,s] = WorkspaceDone(x,v,s)

% ------------------------------------------------------
