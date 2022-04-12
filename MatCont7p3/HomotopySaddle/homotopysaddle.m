function out = homotopysaddle
%
% 
global HTHomds cds
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

[x,x0,p,up,sp,T,eps1,YS,YU] = rearr(arg);
func = BVP_HTHom(x,x0,p,up,sp,T,eps1,YS,YU);  

%------------------------------------------------------
function varargout = jacobian(varargin)
global HTHomds BigJac

[x,x0,p,up,sp,T,eps1,YS,YU] = rearr(varargin{1});

varargout{1} = BVP_HTHom_jac(HTHomds.func,x,x0,p,up,sp,T,eps1,YS,YU);
BigJac = varargout{1};

%-----------------------------------------------------

function hessians(varargin)

%------------------------------------------------------

function varargout = defaultprocessor(varargin)
global HTHomds cds

[x,x0,p,up,sp,T,eps1,YS,YU] = rearr(varargin{1});
v = rearr(varargin{2});

HTHomds.ndim = length(varargin{1});

% update
if ~isempty(HTHomds.ups)
    HTHomds.upold = HTHomds.ups;
end
HTHomds.ups = reshape(x,HTHomds.nphase,HTHomds.tps);
HTHomds.vps = reshape(v,HTHomds.nphase,HTHomds.tps);
if HTHomds.index == 1
    p = HTHomds.P0;  
end
p1 = num2cell(p);
if HTHomds.index == 3
    HTHomds.T = T;
end
for i=1:HTHomds.tps      
    HTHomds.upoldp(:,i) = 2*HTHomds.T*feval(HTHomds.func, 0, HTHomds.ups(:,i),p1{:});
end
HTHomds.eps1 = eps1;
if HTHomds.index == 2 || HTHomds.index == 3
    HTHomds.YS = YS;
    HTHomds.YU = YU;
    HTHomds.x0 = x0;
end

if HTHomds.index == 1
    HTHomds.SParams(HTHomds.ActiveSParams) = sp(HTHomds.ActiveSParams);
elseif HTHomds.index == 2
    HTHomds.SParams(HTHomds.ActiveSParams) = sp;
end
if HTHomds.index == 1
    HTHomds.UParams(HTHomds.ActiveUParams) = up(HTHomds.ActiveUParams);
end

HTHomds.TestTolerance = cds.options.TestTolerance;

% Update dimensions
% -----------------
A = cjac(HTHomds.func,HTHomds.Jacobian,HTHomds.x0,p1,[]);
D = eig(A);

% nneg = dimension of stable subspace
%HTHomds.nneg = sum(real(D) < 0);
% If one eigenvalue is (practically) zero, and the one of the subspaces has
% zero dimension, change this dimension with 1.
%if (HTHomds.nneg == HTHomds.nphase)
%    if min(abs(real(D))) < 1e-9
%        HTHomds.nneg = HTHomds.nneg -1;
%    end
%end
%if (HTHomds.nneg == 0)
%    if min(abs(real(D))) < 1e-9
%        HTHomds.nneg = HTHomds.nneg +1;
%    end
%end
%HTHomds.npos = HTHomds.nphase-HTHomds.nneg;

if nargin > 2
    % set data in special point structure
    s = varargin{3};
    s.data.timemesh = HTHomds.msh;
    s.data.ntst = HTHomds.ntst;
    s.data.ncol = HTHomds.ncol;
    s.data.parametervalues = p1;
    if HTHomds.index == 3
        s.data.T = T;
    end
    varargout{3} = s;
end

% all done succesfully
varargout{1} = 0;
varargout{2} = HTHomds.msh';
  
if (cds.options.Eigenvalues==1)
    varargout{2} = [varargout{2}; D];
end



%-------------------------------------------------------
function option = options
global HTHomds 
  
% Check for symbolic derivatives in odefile
symord = 0; 
if ~isempty(HTHomds.Jacobian), symord = 1; end
if ~isempty(HTHomds.Hessians), symord = 2; end
if ~isempty(HTHomds.Der3), symord = 3; end
if  ~isempty(HTHomds.Der4), symord = 4; end
if ~isempty(HTHomds.Der5), symord = 5; end

option = contset;
option = contset(option, 'SymDerivative', symord);
option = contset(option, 'Workspace', 1);
option = contset(option, 'Locators', zeros(1,13));
symordp = 0;
if  ~isempty(HTHomds.JacobianP), symordp = 1; end
if ~isempty(HTHomds.HessiansP),  symordp = 2; end
option = contset(option, 'SymDerivativeP', symordp);

%------------------------------------------------------  
function [out, failed] = testf(id, x0, v)
global HTHomds

[x,x0,p,up,sp,T,eps1,YS,YU] = rearr(x0);
   
if HTHomds.index == 1 
    sparam = sp(HTHomds.ActiveSParams);%dat zijn de tau's die nul kunnen worden
elseif HTHomds.index == 2    
    sparam = sp;
end
if HTHomds.index == 1 || HTHomds.index == 2
    res = 1;
    for i = 1:length(sparam)
        res = res*sparam(i);
    end
    out = res;    
else
    out = eps1-HTHomds.eps1tol;
end
failed = [];


%-------------------------------------------------------------
function [out, failed] = userf( userinf, id, x, v)
global HTHomds

dim =size(id,2);
failed = [];

[x,x0,p,up,sp,T,eps1,YS,YU] = rearr(x);
if HTHomds.index == 1    
    p = HTHomds.P0;
    x0 = HTHomds.x0;
end

p1 = num2cell(p);
out(dim) = 0;
for i=1:dim
    lastwarn('');    
    if (userinf(i).state==1)    
        out(i)=feval(HTHomds.user{id(i)},0,x0,p1{:});
    else
        out(i)=0;
    end
    if ~isempty(lastwarn)    
        failed = [failed i];
    end
end

%-----------------------------------------------------------------
function [failed,s] = process(id, x, v, s)
global HTHomds

[x,x0,p,up,sp,T,eps1,YS,YU] = rearr(x);

if HTHomds.index == 1 || HTHomds.index == 2
    fprintf('SParam equal to zero \n');
    s.msg  = sprintf('SParam equal to zero');
else
    fprintf('eps1 small enough \n');
    s.msg  = sprintf('eps1 small enough'); 
end
s.label = 'HTHom'; %dat zorgt ervoor dat je als je de testfunctie aanduidt, dat je dan in Connec_Connec terecht komt

failed = 0;   
   
%-------------------------------------------------------------  
function [S,L] = singmat
global HTHomds

if HTHomds.index == 1 || HTHomds.index == 2  
    S = 0;
    L = ['Prod_asp'];
else
    S = 0;
    L = ['disteps1'];
end

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
global cds HTHomds

cds.adapted = 1;

[x,v] = HTHom_adapt_mesh(x,v);

if HTHomds.index == 2 || HTHomds.index == 3
    YU = HTHomds.YU;
    YS = HTHomds.YS;
    
    Q0S = HTHomds.oldStableQ;
    QbS1 = Q0S(:,1:HTHomds.nneg);
    QbS2 = Q0S(:,HTHomds.nneg+1:end);
    Q0U = HTHomds.oldUnstableQ;
    QbU1 = Q0U(:,1:HTHomds.npos);
    QbU2 = Q0U(:,HTHomds.npos+1:end);
    
    [Q1S,S1,R1] = svd(QbS1 + QbS2*YS);
    [Q1U,S1,R1] = svd(QbU1 + QbU2*YU);   
    
    %stable case
    Q1S1=Q1S(:,1:HTHomds.nneg);
    Q1S2=Q1S(:,HTHomds.nneg+1:end);

    Q1S11=Q1S1*((Q1S1'*Q1S1)\(Q1S1'*QbS1));
    for j=1:HTHomds.nneg
        Q1S11(:,j)=Q1S11(:,j)/norm(Q1S11(:,j));
        for i=j+1:HTHomds.nneg
            Q1S11(:,i)=Q1S11(:,i)-Q1S11(:,j)*(Q1S11(:,j)'*Q1S11(:,i));
        end
    end
    Q1S21=Q1S2*((Q1S2'*Q1S2)\(Q1S2'*QbS2));
    for j=1:HTHomds.npos
        Q1S21(:,j)=Q1S21(:,j)/norm(Q1S21(:,j));
        for i=j+1:HTHomds.npos
            Q1S21(:,i)=Q1S21(:,i)-Q1S21(:,j)*(Q1S21(:,j)'*Q1S21(:,i));
        end
    end
    Q1S=[Q1S11 Q1S21];

    %unstable case
    Q1U1=Q1U(:,1:HTHomds.npos);
    Q1U2=Q1U(:,HTHomds.npos+1:end);

    Q1U11=Q1U1*((Q1U1'*Q1U1)\(Q1U1'*QbU1));
    for j=1:HTHomds.npos
        Q1U11(:,j)=Q1U11(:,j)/norm(Q1U11(:,j));
        for i=j+1:HTHomds.npos
            Q1U11(:,i)=Q1U11(:,i)-Q1U11(:,j)*(Q1U11(:,j)'*Q1U11(:,i));
        end
    end
    Q1U21=Q1U2*((Q1U2'*Q1U2)\(Q1U2'*QbU2));
    for j=1:HTHomds.nneg
        Q1U21(:,j)=Q1U21(:,j)/norm(Q1U21(:,j));
        for i=j+1:HTHomds.nneg
            Q1U21(:,i)=Q1U21(:,i)-Q1U21(:,j)*(Q1U21(:,j)'*Q1U21(:,i));
        end
    end
    Q1U=[Q1U11 Q1U21]; 
    
    HTHomds.oldStableQ = Q1S;
    HTHomds.oldUnstableQ = Q1U;

    HTHomds.YS = zeros(size(HTHomds.YS));
    HTHomds.YU = zeros(size(HTHomds.YU));

    x(end-size(HTHomds.YS,1)*size(HTHomds.YS,2)-size(HTHomds.YU,1)*size(HTHomds.YU,2)+1:end)=zeros(size(HTHomds.YS,1)*size(HTHomds.YS,2)+size(HTHomds.YU,1)*size(HTHomds.YU,2),1);
end

res = 1;
% ---------------------------------------------------------------

%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------


%--------------------------------------------------------------
function [x,x0,p,up,sp,T,eps1,YS,YU] = rearr(x1)
% Rearranges x1 into all of its components
global HTHomds

if HTHomds.index == 1
    x = x1(HTHomds.coords);

    up = HTHomds.UParams;
    up(HTHomds.ActiveUParams) = x1(HTHomds.ncoords+1:HTHomds.ncoords+size(HTHomds.ActiveUParams,2));

    sp = HTHomds.SParams;
    sp(HTHomds.ActiveSParams) = x1(HTHomds.ncoords+size(HTHomds.ActiveUParams,2)+1:HTHomds.ncoords+size(HTHomds.ActiveUParams,2)+size(HTHomds.ActiveSParams,2));

    eps1 = x1(HTHomds.ncoords+size(HTHomds.ActiveUParams,2)+size(HTHomds.ActiveSParams,2)+1);
    
    x0 = [];
    p = [];
    T = [];
    YS = [];
    YU = [];
elseif HTHomds.index == 2
    x = x1(HTHomds.coords);    
    x0 = x1(HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.nphase);

    p = HTHomds.P0;
    p(HTHomds.ActiveParams) = x1(HTHomds.PeriodIdx+1);%maar 1 actieve parameter

    sp = x1(HTHomds.PeriodIdx+2);

    eps1 = x1(HTHomds.PeriodIdx+3);

    idx = HTHomds.PeriodIdx+3;

    YU = reshape(x1(idx+1:idx+HTHomds.npos*HTHomds.nneg),HTHomds.nneg,HTHomds.npos);
    idx = idx + HTHomds.npos*HTHomds.nneg;
    YS = reshape(x1(idx+1:idx+HTHomds.npos*HTHomds.nneg),HTHomds.npos,HTHomds.nneg);
    
    up = [];
    T = [];
else 
    x = x1(HTHomds.coords);
    x0 = x1(HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.nphase);

    p = HTHomds.P0;
    p(HTHomds.ActiveParams) = x1(HTHomds.PeriodIdx+1);%1 actieve parameter

    T = x1(HTHomds.PeriodIdx+2);
    eps1 = x1(HTHomds.PeriodIdx+3);

    idx = HTHomds.PeriodIdx+4;

    YU = reshape(x1(idx:idx+HTHomds.npos*HTHomds.nneg-1),HTHomds.nneg,HTHomds.npos);
    idx = idx + HTHomds.npos*HTHomds.nneg;
    YS = reshape(x1(idx:idx+HTHomds.npos*HTHomds.nneg-1),HTHomds.npos,HTHomds.nneg);
    
    up = [];
    sp = [];
end
% -------------------------------------------------------------



% ---------------------------------------------------------------
function WorkspaceInit(x,v)
global HTHomds

HTHomds.cols_p1 = 1:(HTHomds.ncol+1);
HTHomds.cols_p1_coords = 1:(HTHomds.ncol+1)*HTHomds.nphase;
HTHomds.ncol_coord = HTHomds.ncol*HTHomds.nphase;
HTHomds.col_coords = 1:HTHomds.ncol*HTHomds.nphase;
HTHomds.coords = 1:HTHomds.ncoords;
HTHomds.pars = HTHomds.ncoords+1;
HTHomds.tsts = 1:HTHomds.ntst;
HTHomds.cols = 1:HTHomds.ncol;
HTHomds.phases = 1:HTHomds.nphase;
HTHomds.ntstcol = HTHomds.ntst*HTHomds.ncol;

HTHomds.idxmat = reshape(fix((1:((HTHomds.ncol+1)*HTHomds.ntst))/(1+1/HTHomds.ncol))+1,HTHomds.ncol+1,HTHomds.ntst);
HTHomds.dt = HTHomds.msh(HTHomds.tsts+1)-HTHomds.msh(HTHomds.tsts);

HTHomds.wp = kron(HTHomds.wpvec',eye(HTHomds.nphase));
HTHomds.pwwt = kron(HTHomds.wt',eye(HTHomds.nphase));
HTHomds.pwi = HTHomds.wi(ones(1,HTHomds.nphase),:);

HTHomds.wi = nc_weight(HTHomds.ncol)';

[HTHomds.bialt_M1,HTHomds.bialt_M2,HTHomds.bialt_M3,HTHomds.bialt_M4]=bialtaa(HTHomds.nphase);


% ------------------------------------------------------

function [x,v,s] = WorkspaceDone(x,v,s)

% ------------------------------------------------------
