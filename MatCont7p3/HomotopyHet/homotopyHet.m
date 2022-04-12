function out = homotopyHet
%
% 
global HTHetds cds
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

[x,x0,x1,p,up,sp,T,eps1,YS,YU] = rearr(arg);
func = BVP_HTHet(x,x0,x1,p,up,sp,T,eps1,YS,YU);  


%------------------------------------------------------
function varargout = jacobian(varargin)
global HTHetds BigJac

[x,x0,x1,p,up,sp,T,eps1,YS,YU] = rearr(varargin{1});
varargout{1} = BVP_HTHet_jac(HTHetds.func,x,x0,x1,p,up,sp,T,eps1,YS,YU);
BigJac = varargout{1};


%-----------------------------------------------------

function hessians(varargin)

%------------------------------------------------------

function varargout = defaultprocessor(varargin)
global HTHetds cds

[x,x0,x1,p,up,sp,T,eps1,YS,YU] = rearr(varargin{1});
v = rearr(varargin{2});

HTHetds.ndim = length(varargin{1});

% update
if ~isempty(HTHetds.ups)
    HTHetds.upold = HTHetds.ups;
end
HTHetds.ups = reshape(x,HTHetds.nphase,HTHetds.tps);
HTHetds.vps = reshape(v,HTHetds.nphase,HTHetds.tps);
if HTHetds.index == 1
    p = HTHetds.P0;  
end
p1 = num2cell(p);
if HTHetds.index == 3
    HTHetds.T = T;
end
for i=1:HTHetds.tps      
    HTHetds.upoldp(:,i) = 2*HTHetds.T*feval(HTHetds.func, 0, HTHetds.ups(:,i),p1{:});
end

HTHetds.eps1 = eps1;
if HTHetds.index == 2 || HTHetds.index == 3
    HTHetds.YS = YS;
    HTHetds.YU = YU;
    HTHetds.x0 = x0;
    HTHetds.x1 = x1;
end

if HTHetds.index == 1 || HTHetds.index == 2
    HTHetds.SParams(HTHetds.ActiveSParams) = sp(HTHetds.ActiveSParams);
end

if HTHetds.index == 1
    HTHetds.UParams(HTHetds.ActiveUParams) = up(HTHetds.ActiveUParams);
end

HTHetds.TestTolerance = cds.options.TestTolerance;

% Update dimensions
% -----------------
A1 = cjac(HTHetds.func,HTHetds.Jacobian,HTHetds.x0,p1,HTHetds.ActiveParams);
D1 = eig(A1);
A2 = cjac(HTHetds.func,HTHetds.Jacobian,HTHetds.x1,p1,HTHetds.ActiveParams);
D2 = eig(A2);

if nargin > 2
    % set data in special point structure
    s = varargin{3};
    s.data.timemesh = HTHetds.msh;
    s.data.ntst = HTHetds.ntst;
    s.data.ncol = HTHetds.ncol;
    s.data.parametervalues = p1;
    if HTHetds.index == 3
        s.data.T = T;
    end
    varargout{3} = s;
end

% all done succesfully
varargout{1} = 0;
varargout{2} = HTHetds.msh';
  
if (cds.options.Eigenvalues==1)
    varargout{2} = [varargout{2}; D1; D2];
end


%-------------------------------------------------------
function option = options
global HTHetds 

% Check for symbolic derivatives in odefile
symord = 0; 
if ~isempty(HTHetds.Jacobian), symord = 1; end
if ~isempty(HTHetds.Hessians), symord = 2; end
if ~isempty(HTHetds.Der3), symord = 3; end
if  ~isempty(HTHetds.Der4), symord = 4; end
if ~isempty(HTHetds.Der5), symord = 5; end

option = contset;
option = contset(option, 'SymDerivative', symord);
option = contset(option, 'Workspace', 1);
option = contset(option, 'Locators', zeros(1,13));
symordp = 0;
if  ~isempty(HTHetds.JacobianP), symordp = 1; end
if ~isempty(HTHetds.HessiansP),  symordp = 2; end
option = contset(option, 'SymDerivativeP', symordp);

%------------------------------------------------------  
function [out, failed] = testf(id, x0, v)
global HTHetds

[x,x0,x1,p,up,sp,T,eps1,YS,YU] = rearr(x0);
   
if HTHetds.index == 1 || HTHetds.index == 2
    sparam = sp(HTHetds.ActiveSParams);%dat zijn de tau's die nul kunnen worden
end
if HTHetds.index == 1 || HTHetds.index == 2
    res = 1;
    for i = 1:length(sparam)
        res = res*sparam(i);
    end
    out = res;    
else
    out = eps1-HTHetds.eps1tol;
end
failed = [];

%-------------------------------------------------------------
function [out, failed] = userf( userinf, id, x, v)
global HTHetds

dim =size(id,2);
failed = [];
[x,x0,x1,p,up,sp,T,eps1,YS,YU] = rearr(x);
if HTHetds.index == 1    
    p = HTHetds.P0;
    x0 = HTHetds.x0;
end

p1 = num2cell(p);
out(dim) = 0;
for i=1:dim
    lastwarn('');    
    if (userinf(i).state==1)    
        out(i)=feval(HTHetds.user{id(i)},0,x0,p1{:});
    else
        out(i)=0;
    end
    if ~isempty(lastwarn)    
        failed = [failed i];
    end
end

%-----------------------------------------------------------------
function [failed,s] = process(id, x, v, s)
global HTHetds

[x,x0,x1,p,up,sp,T,eps1,YS,YU] = rearr(x);

if HTHetds.index == 1 || HTHetds.index == 2
    fprintf('SParam equal to zero \n');
    s.msg  = sprintf('SParam equal to zero');
else
    fprintf('eps1 small enough \n');
    s.msg  = sprintf('eps1 small enough'); 
end
s.label = 'HTHet';

failed = 0;   
  
%-------------------------------------------------------------  
function [S,L] = singmat
global HTHetds

if HTHetds.index == 1 || HTHetds.index == 2  
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
global cds HTHetds

cds.adapted = 1;

[x,v] = HTHet_adapt_mesh(x,v);

if HTHetds.index == 2 || HTHetds.index == 3
    if ~isempty(HTHetds.ActiveParams)
        YU = HTHetds.YU;
        YS = HTHetds.YS;

        Q0S = HTHetds.oldStableQ;
        QbS1 = Q0S(:,1:HTHetds.nneg);
        QbS2 = Q0S(:,HTHetds.nneg+1:end);
        Q0U = HTHetds.oldUnstableQ;
        QbU1 = Q0U(:,1:HTHetds.npos);
        QbU2 = Q0U(:,HTHetds.npos+1:end);
    
        [Q1S,S1,R1] = svd(QbS1 + QbS2*YS);
        [Q1U,S1,R1] = svd(QbU1 + QbU2*YU);
    
        %stable case
        Q1S1=Q1S(:,1:HTHetds.nneg);
        Q1S2=Q1S(:,HTHetds.nneg+1:end);

        Q1S11=Q1S1*((Q1S1'*Q1S1)\(Q1S1'*QbS1));
        for j=1:HTHetds.nneg
            Q1S11(:,j)=Q1S11(:,j)/norm(Q1S11(:,j));
            for i=j+1:HTHetds.nneg
                Q1S11(:,i)=Q1S11(:,i)-Q1S11(:,j)*(Q1S11(:,j)'*Q1S11(:,i));
            end
        end
        Q1S21=Q1S2*((Q1S2'*Q1S2)\(Q1S2'*QbS2));
        for j=1:HTHetds.nphase-HTHetds.nneg
            Q1S21(:,j)=Q1S21(:,j)/norm(Q1S21(:,j));
            for i=j+1:HTHetds.nphase-HTHetds.nneg
                Q1S21(:,i)=Q1S21(:,i)-Q1S21(:,j)*(Q1S21(:,j)'*Q1S21(:,i));
            end
        end
        Q1S=[Q1S11 Q1S21];

        %unstable case
        Q1U1=Q1U(:,1:HTHetds.npos);
        Q1U2=Q1U(:,HTHetds.npos+1:end);

        Q1U11=Q1U1*((Q1U1'*Q1U1)\(Q1U1'*QbU1));
        for j=1:HTHetds.npos
            Q1U11(:,j)=Q1U11(:,j)/norm(Q1U11(:,j));
            for i=j+1:HTHetds.npos
                Q1U11(:,i)=Q1U11(:,i)-Q1U11(:,j)*(Q1U11(:,j)'*Q1U11(:,i));
            end
        end
        Q1U21=Q1U2*((Q1U2'*Q1U2)\(Q1U2'*QbU2));
        for j=1:HTHetds.nphase-HTHetds.npos
            Q1U21(:,j)=Q1U21(:,j)/norm(Q1U21(:,j));
            for i=j+1:HTHetds.nphase-HTHetds.npos
                Q1U21(:,i)=Q1U21(:,i)-Q1U21(:,j)*(Q1U21(:,j)'*Q1U21(:,i));
            end
        end
        Q1U=[Q1U11 Q1U21];
    
        HTHetds.oldStableQ = Q1S;
        HTHetds.oldUnstableQ = Q1U;

        HTHetds.YS = zeros(size(HTHetds.YS));
        HTHetds.YU = zeros(size(HTHetds.YU));

        x(end-size(HTHetds.YS,1)*size(HTHetds.YS,2)-size(HTHetds.YU,1)*size(HTHetds.YU,2)+1:end)=zeros(size(HTHetds.YS,1)*size(HTHetds.YS,2)+size(HTHetds.YU,1)*size(HTHetds.YU,2),1);
    end
end

res = 1;

% ---------------------------------------------------------------

%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------


%--------------------------------------------------------------
function [x,x0,x1,p,up,sp,T,eps1,YS,YU] = rearr(y)
% Rearranges x1 into all of its components
global HTHetds

if HTHetds.index == 1
    x = y(HTHetds.coords);
    
    up = HTHetds.UParams;
    up(HTHetds.ActiveUParams) = y(HTHetds.ncoords+1:HTHetds.ncoords+size(HTHetds.ActiveUParams,2));

    sp = HTHetds.SParams;
    sp(HTHetds.ActiveSParams) = y(HTHetds.ncoords+size(HTHetds.ActiveUParams,2)+1:HTHetds.ncoords+size(HTHetds.ActiveUParams,2)+size(HTHetds.ActiveSParams,2));

    eps1 = y(HTHetds.ncoords+size(HTHetds.ActiveUParams,2)+size(HTHetds.ActiveSParams,2)+1);
    
    x0 = [];
    x1 = [];
    p = [];
    T = [];
    YS = [];
    YU = [];
  
elseif HTHetds.index == 2
    x = y(HTHetds.coords);    
    x0 = y(HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase);
    x1 = y(HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+2*HTHetds.nphase);

    p = HTHetds.P0;
    p(HTHetds.ActiveParams) = y(HTHetds.PeriodIdx+1:HTHetds.PeriodIdx+length(HTHetds.ActiveParams));
    
    sp = HTHetds.SParams;
    sp(HTHetds.ActiveSParams) = y(HTHetds.PeriodIdx+length(HTHetds.ActiveParams)+1:HTHetds.PeriodIdx+length(HTHetds.ActiveParams)+length(HTHetds.ActiveSParams));

    eps1 = y(HTHetds.PeriodIdx+length(HTHetds.ActiveParams)+length(HTHetds.ActiveSParams)+1);

    idx = HTHetds.PeriodIdx+length(HTHetds.ActiveParams)+length(HTHetds.ActiveSParams)+1;

    YU = reshape(y(idx+1:idx+(HTHetds.nphase-HTHetds.npos)*HTHetds.npos),HTHetds.nphase-HTHetds.npos,HTHetds.npos);
    idx = idx + (HTHetds.nphase-HTHetds.npos)*HTHetds.npos;
    YS = reshape(y(idx+1:idx+(HTHetds.nphase-HTHetds.nneg)*HTHetds.nneg),HTHetds.nphase-HTHetds.nneg,HTHetds.nneg);
    
    up = [];
    T = [];
else 
    x = y(HTHetds.coords);
    x0 = y(HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase);
    x1 = y(HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+2*HTHetds.nphase);

    p = HTHetds.P0;
    p(HTHetds.ActiveParams) = y(HTHetds.PeriodIdx+1:HTHetds.PeriodIdx+length(HTHetds.ActiveParams));

    T = y(HTHetds.PeriodIdx+length(HTHetds.ActiveParams)+1);
    eps1 = y(HTHetds.PeriodIdx+length(HTHetds.ActiveParams)+2);

    idx = HTHetds.PeriodIdx+length(HTHetds.ActiveParams)+3;

    YU = reshape(y(idx:idx+(HTHetds.nphase-HTHetds.npos)*HTHetds.npos-1),HTHetds.nphase-HTHetds.npos,HTHetds.npos);
    idx = idx + (HTHetds.nphase-HTHetds.npos)*HTHetds.npos;
    YS = reshape(y(idx:idx+(HTHetds.nphase-HTHetds.nneg)*HTHetds.nneg-1),HTHetds.nphase-HTHetds.nneg,HTHetds.nneg);
    
    up = [];
    sp = [];
end
% -------------------------------------------------------------



% ---------------------------------------------------------------
function WorkspaceInit(x,v)
global HTHetds

HTHetds.cols_p1 = 1:(HTHetds.ncol+1);
HTHetds.cols_p1_coords = 1:(HTHetds.ncol+1)*HTHetds.nphase;
HTHetds.ncol_coord = HTHetds.ncol*HTHetds.nphase;
HTHetds.col_coords = 1:HTHetds.ncol*HTHetds.nphase;
HTHetds.coords = 1:HTHetds.ncoords;
HTHetds.pars = HTHetds.ncoords+2*HTHetds.nphase+(1:length(HTHetds.ActiveParams));
HTHetds.tsts = 1:HTHetds.ntst;
HTHetds.cols = 1:HTHetds.ncol;
HTHetds.phases = 1:HTHetds.nphase;
HTHetds.ntstcol = HTHetds.ntst*HTHetds.ncol;

HTHetds.idxmat = reshape(fix((1:((HTHetds.ncol+1)*HTHetds.ntst))/(1+1/HTHetds.ncol))+1,HTHetds.ncol+1,HTHetds.ntst);
HTHetds.dt = HTHetds.msh(HTHetds.tsts+1)-HTHetds.msh(HTHetds.tsts);

HTHetds.wp = kron(HTHetds.wpvec',eye(HTHetds.nphase));
HTHetds.pwwt = kron(HTHetds.wt',eye(HTHetds.nphase));
HTHetds.pwi = HTHetds.wi(ones(1,HTHetds.nphase),:);

HTHetds.wi = nc_weight(HTHetds.ncol)';

[HTHetds.bialt_M1,HTHetds.bialt_M2,HTHetds.bialt_M3,HTHetds.bialt_M4]=bialtaa(HTHetds.nphase);

% ------------------------------------------------------

function [x,v,s] = WorkspaceDone(x,v,s)

% ------------------------------------------------------
