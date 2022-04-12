function out = homotopysaddlenode
%
% 
global HTHSNds cds
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

[x,up,sp,T,eps1] = rearr(arg);
func = BVP_HTHSN(x,up,sp,T,eps1);  

%------------------------------------------------------
function varargout = jacobian(varargin)
global HTHSNds BigJac

[x,up,sp,T,eps1] = rearr(varargin{1});

varargout{1} = BVP_HTHSN_jac(HTHSNds.func,x,up,sp,T,eps1);
BigJac = varargout{1};

%-----------------------------------------------------

function hessians(varargin)

%------------------------------------------------------

function varargout = defaultprocessor(varargin)
global HTHSNds cds

[x,up,sp,T,eps1] = rearr(varargin{1});
v = rearr(varargin{2});

HTHSNds.ndim = length(varargin{1});

% update
if ~isempty(HTHSNds.ups)
    HTHSNds.upold = HTHSNds.ups;
end
HTHSNds.ups = reshape(x,HTHSNds.nphase,HTHSNds.tps);
HTHSNds.vps = reshape(v,HTHSNds.nphase,HTHSNds.tps);

p = HTHSNds.P0;  
p1 = num2cell(p);

if HTHSNds.index == 2
    HTHSNds.T = T;
end

for i=1:HTHSNds.tps      
    HTHSNds.upoldp(:,i) = 2*HTHSNds.T*feval(HTHSNds.func, 0, HTHSNds.ups(:,i),p1{:});
end
HTHSNds.eps1 = eps1;

if HTHSNds.index == 1
    HTHSNds.SParams(HTHSNds.ActiveSParams) = sp(HTHSNds.ActiveSParams);
    HTHSNds.UParams(HTHSNds.ActiveUParams) = up(HTHSNds.ActiveUParams);
end

HTHSNds.TestTolerance = cds.options.TestTolerance;

% Update dimensions
% -----------------
A = cjac(HTHSNds.func,HTHSNds.Jacobian,HTHSNds.x0,p1,[]);
D = eig(A);

% nneg = dimension of stable subspace
%HTHSNds.nneg = sum(real(D) < 0);
% If one eigenvalue is (practically) zero, and the one of the subspaces has
% zero dimension, change this dimension with 1.
%if (HTHSNds.nneg == HTHSNds.nphase)
%    if min(abs(real(D))) < 1e-9
%        HTHSNds.nneg = HTHSNds.nneg -1;
%    end
%end
%if (HTHSNds.nneg == 0)
%    if min(abs(real(D))) < 1e-9
%        HTHSNds.nneg = HTHSNds.nneg +1;
%    end
%end
%HTHSNds.npos = HTHSNds.nphase-HTHSNds.nneg;

if nargin > 2
    % set data in special point structure
    s = varargin{3};
    s.data.timemesh = HTHSNds.msh;
    s.data.ntst = HTHSNds.ntst;
    s.data.ncol = HTHSNds.ncol;
    s.data.parametervalues = p1;
    if HTHSNds.index == 2
        s.data.T = T;
    end
    varargout{3} = s;
end

% all done succesfully
varargout{1} = 0;
varargout{2} = HTHSNds.msh';
  
if (cds.options.Eigenvalues==1)
    varargout{2} = [varargout{2}; D];
end

%-------------------------------------------------------
function option = options
global HTHSNds 
  
% Check for symbolic derivatives in odefile
symord = 0; 
if ~isempty(HTHSNds.Jacobian), symord = 1; end
if ~isempty(HTHSNds.Hessians), symord = 2; end
if ~isempty(HTHSNds.Der3), symord = 3; end
if  ~isempty(HTHSNds.Der4), symord = 4; end
if ~isempty(HTHSNds.Der5), symord = 5; end

option = contset;
option = contset(option, 'SymDerivative', symord);
option = contset(option, 'Workspace', 1);
option = contset(option, 'Locators', zeros(1,13));
symordp = 0;
if  ~isempty(HTHSNds.JacobianP), symordp = 1; end
if ~isempty(HTHSNds.HessiansP),  symordp = 2; end
option = contset(option, 'SymDerivativeP', symordp);

%------------------------------------------------------  
function [out, failed] = testf(id, x0, v)
global HTHSNds

[x,up,sp,T,eps1] = rearr(x0);
   
if HTHSNds.index == 1 
    sparam = sp(HTHSNds.ActiveSParams);%dat zijn de tau's die nul kunnen worden
    res = 1;
    for i = 1:length(sparam)
        res = res*sparam(i);
    end
    out = res;    
else
    out = eps1-HTHSNds.eps1tol;
end
failed = [];


%-------------------------------------------------------------
function [out, failed] = userf( userinf, id, x, v)
global HTHSNds

dim =size(id,2);
failed = [];

p = HTHSNds.P0;
x0 = HTHSNds.x0;

p1 = num2cell(p);
out(dim) = 0;
for i=1:dim
    lastwarn('');    
    if (userinf(i).state==1)    
        out(i)=feval(HTHSNds.user{id(i)},0,x0,p1{:});
    else
        out(i)=0;
    end
    if ~isempty(lastwarn)    
        failed = [failed i];
    end
end

%-----------------------------------------------------------------
function [failed,s] = process(id, x, v, s)
global HTHSNds

[x,up,sp,T,eps1] = rearr(x);

if HTHSNds.index == 1
    fprintf('SParam equal to zero');
    s.msg  = sprintf('SParam equal to zero');
else
    fprintf('eps1 small enough');
    s.msg  = sprintf('eps1 small enough'); 
end
s.label = 'HTHSN';

failed = 0;   
   
%-------------------------------------------------------------  
function [S,L] = singmat
global HTHSNds

if HTHSNds.index == 1
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
global cds HTHSNds

cds.adapted = 1;

[x,v] = HTHSN_adapt_mesh(x,v);

% if HTHSNds.index == 2
%     YU = HTHSNds.YU;
%     YS = HTHSNds.YS;
%     
%     if HTHSNds.nneg
%         Q0S = HTHSNds.oldStableQ;
%         QbS1 = Q0S(:,1:HTHSNds.nneg);
%         QbS2 = Q0S(:,HTHSNds.nneg+1:end);
%         [Q1S,S1,R1] = svd(QbS1 + QbS2*YS);
% 
%         Q1S1=Q1S(:,1:HTHSNds.nneg);
%         Q1S2=Q1S(:,HTHSNds.nneg+1:end);
% 
%         Q1S11=Q1S1*((Q1S1'*Q1S1)\(Q1S1'*QbS1));
%         for j=1:HTHSNds.nneg
%             Q1S11(:,j)=Q1S11(:,j)/norm(Q1S11(:,j));
%             for i=j+1:HTHSNds.nneg
%                 Q1S11(:,i)=Q1S11(:,i)-Q1S11(:,j)*(Q1S11(:,j)'*Q1S11(:,i));
%             end
%         end
%         Q1S21=Q1S2*((Q1S2'*Q1S2)\(Q1S2'*QbS2));
%         for j=1:HTHSNds.nphase-HTHSNds.nneg
%             Q1S21(:,j)=Q1S21(:,j)/norm(Q1S21(:,j));
%             for i=j+1:HTHSNds.nphase-HTHSNds.nneg
%                 Q1S21(:,i)=Q1S21(:,i)-Q1S21(:,j)*(Q1S21(:,j)'*Q1S21(:,i));
%             end
%         end
%         Q1S=[Q1S11 Q1S21];
%        
%         HTHSNds.oldStableQ = Q1S;
%     end
%         
%     if HTHSNds.npos
%         Q0U = HTHSNds.oldUnstableQ;
%         QbU1 = Q0U(:,1:HTHSNds.npos);
%         QbU2 = Q0U(:,HTHSNds.npos+1:end);        
%         [Q1U,S1,R1] = svd(QbU1 + QbU2*YU);    
%         
%         Q1U1=Q1U(:,1:HTHSNds.npos);
%         Q1U2=Q1U(:,HTHSNds.npos+1:end);
%     
%         Q1U11=Q1U1*((Q1U1'*Q1U1)\(Q1U1'*QbU1));
%         for j=1:HTHSNds.npos
%             Q1U11(:,j)=Q1U11(:,j)/norm(Q1U11(:,j));
%             for i=j+1:HTHSNds.npos
%                 Q1U11(:,i)=Q1U11(:,i)-Q1U11(:,j)*(Q1U11(:,j)'*Q1U11(:,i));
%             end
%         end
%         Q1U21=Q1U2*((Q1U2'*Q1U2)\(Q1U2'*QbU2));
%         for j=1:HTHSNds.nhase-HTHSNds.npos
%             Q1U21(:,j)=Q1U21(:,j)/norm(Q1U21(:,j));
%             for i=j+1:HTHSNds.nphase-HTHSNds.npos
%                 Q1U21(:,i)=Q1U21(:,i)-Q1U21(:,j)*(Q1U21(:,j)'*Q1U21(:,i));
%             end
%         end
%         Q1U=[Q1U11 Q1U21];
% 
%         HTHSNds.oldUnstableQ = Q1U;
%     end
%     
%     HTHSNds.YS = zeros(size(HTHSNds.YS));
%     HTHSNds.YU = zeros(size(HTHSNds.YU));
% 
%     x(end-size(HTHSNds.YS,1)*size(HTHSNds.YS,2)-size(HTHSNds.YU,1)*size(HTHSNds.YU,2)+1:end)=zeros(size(HTHSNds.YS,1)*size(HTHSNds.YS,2)+size(HTHSNds.YU,1)*size(HTHSNds.YU,2),1);
%     
%     % Update saddle-node borders   
%     jac = cjac(HTHSNds.func,HTHSNds.Jacobian,HTHSNds.x0,num2cell(HTHSNds.P0),[]);
%     Bord = [   jac       HTHSNds.wvector;...
%            HTHSNds.vvector'      0      ];
%     bunit = [zeros(HTHSNds.nphase,1);1];
%     vext = Bord \ bunit;
%     wext = Bord' \ bunit;
%     %ERROR OR WARNING
%     HTHSNds.vvector = vext(1:HTHSNds.nphase)/norm(vext(1:HTHSNds.nphase));
%     HTHSNds.wvector = wext(1:HTHSNds.nphase)/norm(wext(1:HTHSNds.nphase));
% end

res = 1;
% ---------------------------------------------------------------

%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------


%--------------------------------------------------------------
function [x,up,sp,T,eps1] = rearr(x1)
% Rearranges x1 into all of its components
global HTHSNds

if HTHSNds.index == 1
    x = x1(HTHSNds.coords);

    up = HTHSNds.UParams;
    up(HTHSNds.ActiveUParams) = x1(HTHSNds.ncoords+1:HTHSNds.ncoords+size(HTHSNds.ActiveUParams,2));

    sp = HTHSNds.SParams;
    sp(HTHSNds.ActiveSParams) = x1(HTHSNds.ncoords+size(HTHSNds.ActiveUParams,2)+1:HTHSNds.ncoords+size(HTHSNds.ActiveUParams,2)+size(HTHSNds.ActiveSParams,2));

    eps1 = x1(HTHSNds.ncoords+size(HTHSNds.ActiveUParams,2)+size(HTHSNds.ActiveSParams,2)+1);
    
    T = [];
    
else 
    x = x1(HTHSNds.coords);

    T = x1(HTHSNds.ncoords+1);
    eps1 = x1(HTHSNds.ncoords+2);

    up = [];
    sp = [];
end
% -------------------------------------------------------------



% ---------------------------------------------------------------
function WorkspaceInit(x,v)
global HTHSNds

HTHSNds.cols_p1 = 1:(HTHSNds.ncol+1);
HTHSNds.cols_p1_coords = 1:(HTHSNds.ncol+1)*HTHSNds.nphase;
HTHSNds.ncol_coord = HTHSNds.ncol*HTHSNds.nphase;
HTHSNds.col_coords = 1:HTHSNds.ncol*HTHSNds.nphase;
HTHSNds.coords = 1:HTHSNds.ncoords;
HTHSNds.pars = HTHSNds.ncoords+1;
HTHSNds.tsts = 1:HTHSNds.ntst;
HTHSNds.cols = 1:HTHSNds.ncol;
HTHSNds.phases = 1:HTHSNds.nphase;
HTHSNds.ntstcol = HTHSNds.ntst*HTHSNds.ncol;

HTHSNds.idxmat = reshape(fix((1:((HTHSNds.ncol+1)*HTHSNds.ntst))/(1+1/HTHSNds.ncol))+1,HTHSNds.ncol+1,HTHSNds.ntst);
HTHSNds.dt = HTHSNds.msh(HTHSNds.tsts+1)-HTHSNds.msh(HTHSNds.tsts);

HTHSNds.wp = kron(HTHSNds.wpvec',eye(HTHSNds.nphase));
HTHSNds.pwwt = kron(HTHSNds.wt',eye(HTHSNds.nphase));
HTHSNds.pwi = HTHSNds.wi(ones(1,HTHSNds.nphase),:);

HTHSNds.wi = nc_weight(HTHSNds.ncol)';

[HTHSNds.bialt_M1,HTHSNds.bialt_M2,HTHSNds.bialt_M3,HTHSNds.bialt_M4]=bialtaa(HTHSNds.nphase);


% ------------------------------------------------------

function [x,v,s] = WorkspaceDone(x,v,s)

% ------------------------------------------------------
