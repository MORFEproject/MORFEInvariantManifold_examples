function out = homoclinicsaddlenode
%
% homoclinic curve definition file for a problem in odefile
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

  [x,x0,p,T,eps0,eps1,YS,YU] = rearr(arg);
  func = BVP_HSN(x,x0,p,T,eps0,eps1,YS,YU);
  
%------------------------------------------------------

function varargout = jacobian(varargin)

global homds

  [x,x0,p,T,eps0,eps1,YS,YU] = rearr(varargin{1});
  varargout{1} = BVP_HSN_jac(homds.func,x,x0,p,T,eps0,eps1,YS,YU);  
      
  
  
%-----------------------------------------------------

function hessians(varargin)

%------------------------------------------------------

function varargout = defaultprocessor(varargin)
global homds cds

  [x,x0,p,T,eps0,eps1,YS,YU] = rearr(varargin{1});
  if  size(varargin{1}) ==   size(varargin{2})
      v = rearr(varargin{2});
  end
  
  % update
  if ~isempty(homds.ups)
    homds.upold = homds.ups;
  end
  homds.ups = reshape(x,homds.nphase,homds.tps);
  homds.vps = reshape(v,homds.nphase,homds.tps);
  % update upoldp
  p1 = num2cell(p);
  for i=1:homds.tps
    homds.upoldp(:,i) = 2*T*feval(homds.func, 0, homds.ups(:,i), p1{:});
  end
  homds.eps0 = eps0;
  homds.eps1 = eps1;
  homds.YS = YS;
  homds.YU = YU;
  homds.T = T;
  homds.x0 = x0;
  
% Update dimensions
% -----------------
p = num2cell(p);
A = cjac(homds.func,homds.Jacobian,x0,p,homds.ActiveParams);
D = eig(A);
% nneg = dimension of stable subspace
[y,i] = min(abs(D));
homds.nneg = sum(real(D) < 0);
homds.npos = sum(real(D) > 0);
if (homds.nneg+homds.npos == homds.nphase) && (D(i) < 0)
    homds.nneg = homds.nneg-1;
elseif (homds.nneg+homds.npos == homds.nphase)
    homds.npos = homds.npos-1;
end
homds.Ysize = (homds.nneg+1)*homds.npos + homds.nneg*(homds.npos+1);
  
% Normalize bases
% ---------------
  if nargin > 2
    % set data in special point structure
    s = varargin{3};
    s.data.timemesh = homds.msh;
    s.data.ntst = homds.ntst;
    s.data.ncol = homds.ncol;
    s.data.parametervalues = p;
    s.data.T = T;
    varargout{3} = s;
  end
  % all done succesfully
  varargout{1} = 0;
  varargout{2} = homds.msh';
  
  if (cds.options.Eigenvalues==1)
      varargout{2} = [varargout{2}; D];
  end

%-------------------------------------------------------
  
function option = options
global homds 
  % Check for symbolic derivatives in odefile
  
  symord = 0; 
  if ~isempty(homds.Jacobian), symord = 1; end
  if ~isempty(homds.Hessians), symord = 2; end
  if ~isempty(homds.Der3), symord = 3; end
  if ~isempty(homds.Der4), symord = 4; end
  if ~isempty(homds.Der5), symord = 5; end

  option = contset;
  option = contset(option, 'SymDerivative', symord);
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', [0]);
  symordp = 0;
  if ~isempty(homds.JacobianP), symordp = 1; end
  if ~isempty(homds.HessiansP),  symordp = 2; end
  option = contset(option, 'SymDerivativeP', symordp);
  
  
%------------------------------------------------------  
  
function [out, failed] = testf(id, x0, v)
global homds 

[x,x0,p,T,eps0,eps1,YS,YU] = rearr(x0);
ups = reshape(x,homds.nphase,homds.tps);
A = cjac(homds.func,homds.Jacobian,x0,num2cell(p),homds.ActiveParams);

% Non-Central Homoclinic to saddle-node
opts.disp = 0;
[V,D] = eig(A');
[Y,i] = min(abs(diag(D)));
V = V(:,i);

psi1 = 1/homds.T * ((ups(:,1) - x0)' * V);
psi2 = 1/homds.T * ((ups(:,end) - x0)' * V);
res(1) = psi1 * psi2;

out = res;
failed = [];

%-------------------------------------------------------------

function [out, failed] = userf(userinf, id, x, v)
global  homds
dim =size(id,2);
failed = [];
out(dim) = 0;
for i=1:dim
  lastwarn('');
  [x,x0,p,T,eps0,eps1,YS,YU] = rearr(x); p = num2cell(p);
  if (userinf(i).state==1)
      out(i)=feval(homds.user{id(i)},0,x0,p{:});
  else
      out(i)=0;
  end
  if ~isempty(lastwarn)
    failed = [failed i];
  end
end

%-----------------------------------------------------------------

function [failed,s] = process(id, x, v, s)
global homds
[x,x0,p,T,eps0,eps1,YS,YU] = rearr(x);
AP = p(homds.ActiveParams);
switch id
    case 1
        fprintf('Non-Central Homoclinic to Saddle-Node, parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Non-Central Homoclinic to Saddle-Node');
end
failed = 0;

%-------------------------------------------------------------  

function [S,L] = singmat

  S =  0 ;
  L = 'NCH' ;

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
WorkspaceDone

%-----------------------------------------------------------

function [res,x,v] = adapt(x,v)
global homds

cds.adapted = 1;

YU = homds.YU;
YS = homds.YS;

[x,v] = Hom_adapt_mesh(x,v);

if homds.nneg
    Q0S = homds.oldStableQ;    
    QbS1 = Q0S(:,1:homds.nneg);
    QbS2 = Q0S(:,homds.nneg+1:end);    
    [Q1S,S1,R1] = svd(QbS1 + QbS2*YS);

    Q1S1=Q1S(:,1:homds.nneg);
    Q1S2=Q1S(:,homds.nneg+1:end);

    Q1S11=Q1S1*((Q1S1'*Q1S1)\(Q1S1'*QbS1));
    for j=1:homds.nneg    
        Q1S11(:,j)=Q1S11(:,j)/norm(Q1S11(:,j));
        for i=j+1:homds.nneg
            Q1S11(:,i)=Q1S11(:,i)-Q1S11(:,j)*(Q1S11(:,j)'*Q1S11(:,i));
        end
    end
    Q1S21=Q1S2*((Q1S2'*Q1S2)\(Q1S2'*QbS2));
    for j=1:homds.nphase-homds.nneg
        Q1S21(:,j)=Q1S21(:,j)/norm(Q1S21(:,j));
        for i=j+1:homds.nphase-homds.nneg
            Q1S21(:,i)=Q1S21(:,i)-Q1S21(:,j)*(Q1S21(:,j)'*Q1S21(:,i));
        end
    end
    Q1S=[Q1S11 Q1S21];

    homds.oldStableQ = Q1S;
end

if homds.npos
    Q0U = homds.oldUnstableQ;
    QbU1 = Q0U(:,1:homds.npos);
    QbU2 = Q0U(:,homds.npos+1:end);
    [Q1U,S1,R1] = svd(QbU1 + QbU2*YU);

    Q1U1=Q1U(:,1:homds.npos);
    Q1U2=Q1U(:,homds.npos+1:end);

    Q1U11=Q1U1*((Q1U1'*Q1U1)\(Q1U1'*QbU1));
    for j=1:homds.npos
        Q1U11(:,j)=Q1U11(:,j)/norm(Q1U11(:,j));
        for i=j+1:homds.npos
            Q1U11(:,i)=Q1U11(:,i)-Q1U11(:,j)*(Q1U11(:,j)'*Q1U11(:,i));
        end
    end
    Q1U21=Q1U2*((Q1U2'*Q1U2)\(Q1U2'*QbU2));
    for j=1:homds.nphase-homds.npos
        Q1U21(:,j)=Q1U21(:,j)/norm(Q1U21(:,j));
        for i=j+1:homds.nphase-homds.npos
            Q1U21(:,i)=Q1U21(:,i)-Q1U21(:,j)*(Q1U21(:,j)'*Q1U21(:,i));
        end
    end
    Q1U=[Q1U11 Q1U21];

    homds.oldUnstableQ = Q1U;
end

homds.YS = zeros(size(homds.YS));
homds.YU = zeros(size(homds.YU));

x(end-size(homds.YS,1)*size(homds.YS,2)-size(homds.YU,1)*size(homds.YU,2)+1:end)=zeros(size(homds.YS,1)*size(homds.YS,2)+size(homds.YU,1)*size(homds.YU,2),1);

% Update saddle-node borders
[xtmp,x0,p,T,eps0,eps1,YS,YU] = rearr(x);
jac = cjac(homds.func,homds.Jacobian,x0,num2cell(p),homds.ActiveParams);
Bord = [   jac       homds.wvector;...
       homds.vvector'      0      ];
bunit = [zeros(homds.nphase,1);1];
vext = Bord \ bunit;
wext = Bord' \ bunit;
%ERROR OR WARNING
homds.vvector = vext(1:homds.nphase)/norm(vext(1:homds.nphase));
homds.wvector = wext(1:homds.nphase)/norm(wext(1:homds.nphase));

res = 1;

%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------
 
function [x,x0,p,T,eps0,eps1,YS,YU] = rearr(x1)
% Rearranges x1 into all of its components
global homds

x = x1(homds.coords);
x0 = x1(homds.ncoords+1:homds.ncoords+homds.nphase);

p = homds.P0;
% p(homds.ActiveParams) = x1(homds.ncoords+homds.nphase+1:homds.ncoords+homds.nphase+2);
% idx = homds.ncoords+homds.nphase+3;
p(homds.ActiveParams) = x1(homds.PeriodIdx+1:homds.PeriodIdx+2);
idx = homds.PeriodIdx+3;

if homds.extravec(1)
    T = x1(idx);
    idx = idx+1;
else
    T = homds.T;
end
if homds.extravec(2)
    eps0 = x1(idx);
    idx = idx+1;
else
    eps0 = homds.eps0;
end
if homds.extravec(3)
    eps1 = x1(idx);
    idx = idx+1;
else
    eps1 = homds.eps1;
end
    if homds.npos
        if ~homds.nneg
            YU = reshape(x1(end-(homds.nneg+1)*homds.npos+1:end),homds.nneg+1,homds.npos);
        else
            YU = reshape(x1(end-(homds.nneg+1)*homds.npos-(homds.npos+1)*homds.nneg+1:end-(homds.npos+1)*homds.nneg),homds.nneg+1,homds.npos);
        end
        idx = idx + homds.npos*homds.nneg;
    else
        YU = [];
    end
    if homds.nneg
        YS = reshape(x1(end-(homds.npos+1)*homds.nneg+1:end),homds.npos+1,homds.nneg);
    else
        YS = [];
    end
    
% ---------------------------------------------------------------

function WorkspaceInit(x,v)
global homds
homds.cols_p1 = 1:(homds.ncol+1);
homds.cols_p1_coords = 1:(homds.ncol+1)*homds.nphase;
homds.ncol_coord = homds.ncol*homds.nphase;
homds.col_coords = 1:homds.ncol*homds.nphase;
homds.coords = 1:homds.ncoords;
homds.pars = homds.ncoords+(1:2);
homds.tsts = 1:homds.ntst;
homds.cols = 1:homds.ncol;
homds.phases = 1:homds.nphase;
homds.ntstcol = homds.ntst*homds.ncol;

homds.idxmat = reshape(fix((1:((homds.ncol+1)*homds.ntst))/(1+1/homds.ncol))+1,homds.ncol+1,homds.ntst);
homds.dt = homds.msh(homds.tsts+1)-homds.msh(homds.tsts);

homds.wp = kron(homds.wpvec',eye(homds.nphase));
homds.pwwt = kron(homds.wt',eye(homds.nphase));
homds.pwi = homds.wi(ones(1,homds.nphase),:);

homds.wi = nc_weight(homds.ncol)';

[homds.bialt_M1,homds.bialt_M2,homds.bialt_M3,homds.bialt_M4]=bialtaa(homds.nphase);

% ------------------------------------------------------

function [x,v,s] = WorkspaceDone(x,v,s)

%------------------------------------------------------------
