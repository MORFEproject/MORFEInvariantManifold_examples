function [x,v] = init_HTHSN_HTHSN(odefile, x, v, s, p, up, aup, sp, asp, ntst, ncol,T,eps1,eps1tol)

global HTHSNds cds

oldHTHSNds = HTHSNds;
% initialize HTHSNds

[x1,v] = init_HTHSNds(odefile,x,v,p,up,aup,sp,asp,ntst,ncol,T,eps1,eps1tol,s);
func_handles = feval(odefile);
symord = 0; 
symordp = 0;

if     ~isempty(func_handles{9}),   symord = 5; 
elseif ~isempty(func_handles{8}),   symord = 4; 
elseif ~isempty(func_handles{7}),   symord = 3; 
elseif ~isempty(func_handles{5}),   symord = 2; 
elseif ~isempty(func_handles{3}),   symord = 1; 
end
if     ~isempty(func_handles{6}),   symordp = 2; 
elseif ~isempty(func_handles{4}),   symordp = 1; 
end
if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.symjac = 1;
cds.symhess = 0;


HTHSNds.odefile = odefile;
HTHSNds.func = func_handles{2};
HTHSNds.Jacobian  = func_handles{3};
HTHSNds.JacobianP = func_handles{4};
HTHSNds.Hessians  = func_handles{5};
HTHSNds.HessiansP = func_handles{6};
HTHSNds.Der3 = func_handles{7};
HTHSNds.Der4 = func_handles{8};
HTHSNds.Der5 = func_handles{9};


A = cjac(HTHSNds.func,HTHSNds.Jacobian,HTHSNds.x0,num2cell(p),[]);
D = eig(A);

% nneg = dimension of stable subspace
[y,i] = min(abs(D));
B = D([1:i-1 i+1:end]);
HTHSNds.nneg = sum(real(B) < 0);
if (HTHSNds.nneg == HTHSNds.nphase-1)
    if min(abs(real(B))) < 1e-10
        HTHSNds.nneg = HTHSNds.nneg -1;        
    end
end
if (HTHSNds.nneg == 0)
    if min(abs(real(B))) < 1e-10
        HTHSNds.nneg = HTHSNds.nneg +1;
    end
end
HTHSNds.npos = HTHSNds.nphase-HTHSNds.nneg-1;


% COMPOSE X0
% % ----------
if HTHSNds.index == 1
    
    % 1. cycle 
    ups = reshape(x1,HTHSNds.nphase,HTHSNds.tps);
    HTHSNds.upold = ups;
    HTHSNds.PeriodIdx = length(x1);
    % 2. free UParams
    x1 = [x1; HTHSNds.UParams(HTHSNds.ActiveUParams)];
    % 3. free SParams
    x1 = [x1; HTHSNds.SParams(HTHSNds.ActiveSParams)];
    % 4. eps1
    x1 = [x1; HTHSNds.eps1];
    
else
    
    % 1. cycle 
    ups = reshape(x1,HTHSNds.nphase,HTHSNds.tps);
    HTHSNds.upold = ups;
    
    HTHSNds.PeriodIdx = length(x1);
    % 4. T
    x1 = [x1; HTHSNds.T];
    % 5. eps1
    x1 = [x1; HTHSNds.eps1];    
end

x = x1;
v = [];


% ASSIGN SOME VALUES TO HOMOCLINIC FIELDS
% ---------------------------------------

HTHSNds.YS = zeros(HTHSNds.nphase-HTHSNds.nneg,HTHSNds.nneg);
HTHSNds.YU = zeros(HTHSNds.nphase-HTHSNds.npos,HTHSNds.npos);

% Third parameter = unstable_flag, 
% 1 if we want the unstable space, 0 if we want the stable one
[QS, eigvlS] = computeBaseHTHSN(A,1,HTHSNds.nneg);
[QU, eigvlU] = computeBaseHTHSN(A,0,HTHSNds.npos);

HTHSNds.oldStableQ = QS;
HTHSNds.oldUnstableQ = QU;

HTHSNds.ups = [];
HTHSNds.upold = [];
HTHSNds.upoldp = [];



%-----------------------------------------------------------------
function [x1,v1] = init_HTHSNds(odefile,x,v,p,up,aup,sp,asp,ntst,ncol,T,eps1,eps1tol,s)
global HTHSNds

oldHTHSNds = HTHSNds;
HTHSNds = [];

HTHSNds.odefile = odefile;
func_handles = feval(HTHSNds.odefile);
HTHSNds.func = func_handles{2};
HTHSNds.Jacobian  = func_handles{3};
HTHSNds.JacobianP = func_handles{4};
HTHSNds.Hessians  = func_handles{5};
HTHSNds.HessiansP = func_handles{6};
HTHSNds.Der3=[];
siz = size(func_handles,2);
if siz > 9
    j=1;
    for k=10:siz
        HTHSNds.user{j}= func_handles{k};
        j=j+1;
    end
else
    HTHSNds.user=[];
end

HTHSNds.nphase = oldHTHSNds.nphase;
HTHSNds.x0 = x(oldHTHSNds.nphase*(oldHTHSNds.ntst*oldHTHSNds.ncol+1)+1:oldHTHSNds.nphase*(oldHTHSNds.ntst*oldHTHSNds.ncol+1)+oldHTHSNds.nphase);
HTHSNds.ActiveUParams = aup;
HTHSNds.ActiveSParams = asp;
HTHSNds.P0 = p;
HTHSNds.UParams = up;
HTHSNds.SParams = sp;
HTHSNds.index = oldHTHSNds.index;

%%%
HTHSNds.ntst = ntst;
HTHSNds.ncol = ncol;
HTHSNds.tps = ntst*ncol+1;
HTHSNds.ncoords = HTHSNds.tps*oldHTHSNds.nphase;
HTHSNds.coords = 1:HTHSNds.ncoords;
HTHSNds.msh = s.data.timemesh;
x1 = x(1:(oldHTHSNds.ntst*oldHTHSNds.ncol+1)*oldHTHSNds.nphase);

if isempty(v)
    v1 = [];
else
    v1 = v(1:(oldHTHSNds.ntst*oldHTHSNds.ncol+1)*oldHTHSNds.nphase);
end
[x1,v1] = HTHSN_new_mesh(x1,v1,oldHTHSNds.ntst,oldHTHSNds.ncol);

HTHSNds.T = T;
HTHSNds.eps0 = oldHTHSNds.eps0;
HTHSNds.eps1 = eps1;

HTHSNds.TestTolerance = oldHTHSNds.TestTolerance;
HTHSNds.eps1tol = eps1tol;

HTHSNds.cols_p1 = 1:(HTHSNds.ncol+1);
HTHSNds.cols_p1_coords = 1:(HTHSNds.ncol+1)*HTHSNds.nphase;
HTHSNds.ncol_coord = HTHSNds.ncol*HTHSNds.nphase;
HTHSNds.col_coords = 1:HTHSNds.ncol*HTHSNds.nphase;
HTHSNds.pars = HTHSNds.ncoords+(1:3);
HTHSNds.phases = 1:HTHSNds.nphase;
HTHSNds.ntstcol = HTHSNds.ntst*HTHSNds.ncol;
HTHSNds.wp = kron(HTHSNds.wpvec',eye(HTHSNds.nphase));
HTHSNds.pwwt = kron(HTHSNds.wt',eye(HTHSNds.nphase));
HTHSNds.pwi = HTHSNds.wi(ones(1,HTHSNds.nphase),:);

HTHSNds.bialt_M1 = [];
HTHSNds.bialt_M2 = [];
HTHSNds.bialt_M3 = [];
HTHSNds.bialt_M4 = [];
HTHSNds.multipliers = nan;
HTHSNds.monodromy = [];
HTHSNds.multi_r1 = [];
HTHSNds.multi_r2 = [];
HTHSNds.ups = [];
HTHSNds.vps = [];
HTHSNds.tsts = 1:HTHSNds.ntst;
HTHSNds.cols = 1:HTHSNds.ncol;
