function [x,v] = init_HTHom_HTHom(odefile, x, v, s, p, ap, up, aup, sp, asp, ntst, ncol,T,eps1,eps1tol)

global HTHomds cds

oldHTHomds = HTHomds;
% initialize HTHomds

[x1,v] = init_HTHomds(odefile,x,v,p,ap,up,aup,sp,asp,ntst,ncol,T,eps1,eps1tol,s);
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


HTHomds.odefile = odefile;
HTHomds.func = func_handles{2};
HTHomds.Jacobian  = func_handles{3};
HTHomds.JacobianP = func_handles{4};
HTHomds.Hessians  = func_handles{5};
HTHomds.HessiansP = func_handles{6};
HTHomds.Der3 = func_handles{7};
HTHomds.Der4 = func_handles{8};
HTHomds.Der5 = func_handles{9};


A = cjac(HTHomds.func,HTHomds.Jacobian,HTHomds.x0,num2cell(p),HTHomds.ActiveParams);
D = eig(A);

% nneg = dimension of stable subspace
HTHomds.nneg = sum(real(D) < 0);
if (HTHomds.nneg == HTHomds.nphase)
    if min(abs(real(D))) < 1e-10
        HTHomds.nneg = HTHomds.nneg -1;
    end
end
if (HTHomds.nneg == 0)
    if min(abs(real(D))) < 1e-10
        HTHomds.nneg = HTHomds.nneg +1;
    end
end
HTHomds.npos = HTHomds.nphase-HTHomds.nneg;
HTHomds.Ysize = HTHomds.nneg*HTHomds.npos;



% COMPOSE X0
% % ----------
if HTHomds.index == 1    
    % 1. cycle 
    ups = reshape(x1,HTHomds.nphase,HTHomds.tps);
    HTHomds.upold = ups;
    HTHomds.PeriodIdx = length(x1);
    % 2. free UParams
    x1 = [x1; HTHomds.UParams(HTHomds.ActiveUParams)];
    % 3. free SParams
    x1 = [x1; HTHomds.SParams(HTHomds.ActiveSParams)];
    % 4. eps1
    x1 = [x1; HTHomds.eps1];
    
elseif HTHomds.index == 2
    
    % 1. cycle 
    ups = reshape(x1,HTHomds.nphase,HTHomds.tps);
    HTHomds.upold = ups;
    % 2. equilibrium coordinates
    x1 = [x1; HTHomds.x0];
    HTHomds.PeriodIdx = length(x1);
    % 3. one free parameters
    x1 = [x1; HTHomds.P0(ap)];
    % 4. one free SParams
    x1 = [x1; HTHomds.SParams(HTHomds.ActiveSParams)];
    % 5. eps1
    x1 = [x1; HTHomds.eps1];
    % 6. YS and YU, initialized to 0
    %eerst YU, dan YS
    for i=1:HTHomds.nneg
        x1 = [x1; zeros(HTHomds.npos,1)];
    end
    for i=1:HTHomds.npos
        x1 = [x1; zeros(HTHomds.nneg,1)];
    end

else
    
    % 1. cycle 
    ups = reshape(x1,HTHomds.nphase,HTHomds.tps);
    HTHomds.upold = ups;
    % 2. equilibrium coordinates
    x1 = [x1; HTHomds.x0];
    HTHomds.PeriodIdx = length(x1);
    % 3. one free parameter
    x1 = [x1; HTHomds.P0(ap)];
    % 4. T
    x1 = [x1; HTHomds.T];
    % 5. eps1
    x1 = [x1; HTHomds.eps1];
    % 6. YS and YU, initialized to 0
    for i=1:HTHomds.nneg
        x1 = [x1; zeros(HTHomds.npos,1)];
    end
    for i=1:HTHomds.npos
        x1 = [x1; zeros(HTHomds.nneg,1)];
    end
    
end

x = x1;
v = [];


% ASSIGN SOME VALUES TO HOMOCLINIC FIELDS
% ---------------------------------------

HTHomds.YS = zeros(HTHomds.npos,HTHomds.nneg);
HTHomds.YU = zeros(HTHomds.nneg,HTHomds.npos);

% Third parameter = unstable_flag, 
% 1 if we want the unstable space, 0 if we want the stable one
[QS, eigvlS] = computeBaseHTHom(A,1,HTHomds.nneg);
[QU, eigvlU] = computeBaseHTHom(A,0,HTHomds.npos);

HTHomds.oldStableQ = QS;
HTHomds.oldUnstableQ = QU;

HTHomds.ups = [];
HTHomds.upold = [];
HTHomds.upoldp = [];


%-----------------------------------------------------------------
function [x1,v1] = init_HTHomds(odefile,x,v,p,ap,up,aup,sp,asp,ntst,ncol,T,eps1,eps1tol,s)
global HTHomds

oldHTHomds = HTHomds;
HTHomds = [];

HTHomds.odefile = odefile;
func_handles = feval(HTHomds.odefile);
HTHomds.func = func_handles{2};
HTHomds.Jacobian  = func_handles{3};
HTHomds.JacobianP = func_handles{4};
HTHomds.Hessians  = func_handles{5};
HTHomds.HessiansP = func_handles{6};
HTHomds.Der3=[];
siz = size(func_handles,2);
if siz > 9
    j=1;
    for k=10:siz
        HTHomds.user{j}= func_handles{k};
        j=j+1;
    end
else
    HTHomds.user=[];
end

HTHomds.nphase = oldHTHomds.nphase;
HTHomds.x0 = x(oldHTHomds.nphase*(oldHTHomds.ntst*oldHTHomds.ncol+1)+1:oldHTHomds.nphase*(oldHTHomds.ntst*oldHTHomds.ncol+1)+oldHTHomds.nphase);

HTHomds.ActiveParams = ap;
HTHomds.ActiveUParams = aup;
HTHomds.ActiveSParams = asp;
HTHomds.P0 = p;
HTHomds.UParams = up;
HTHomds.SParams = sp;
HTHomds.index = oldHTHomds.index;

HTHomds.ntst = ntst;
HTHomds.ncol = ncol;
HTHomds.tps = ntst*ncol+1;
HTHomds.ncoords = HTHomds.tps*HTHomds.nphase;
HTHomds.coords = 1:HTHomds.ncoords;
HTHomds.msh = s.data.timemesh;

x1 = x(1:(oldHTHomds.ntst*oldHTHomds.ncol+1)*oldHTHomds.nphase);
if isempty(v)
    v1 = [];
else
    v1 = v(1:(oldHTHomds.ntst*oldHTHomds.ncol+1)*oldHTHomds.nphase);
end

[x1,v1] = HTHom_new_mesh(x1,v1,oldHTHomds.ntst,oldHTHomds.ncol);
% HTHom_set_ntst_ncol(HTHomds.ntst,HTHomds.ncol,HTHomds.msh);


HTHomds.T = T;
HTHomds.eps0 = oldHTHomds.eps0;
HTHomds.eps1 = eps1;

HTHomds.TestTolerance = oldHTHomds.TestTolerance;
HTHomds.eps1tol = eps1tol;

HTHomds.cols_p1 = 1:(HTHomds.ncol+1);
HTHomds.cols_p1_coords = 1:(HTHomds.ncol+1)*HTHomds.nphase;
HTHomds.ncol_coord = HTHomds.ncol*HTHomds.nphase;
HTHomds.col_coords = 1:HTHomds.ncol*HTHomds.nphase;
HTHomds.pars = HTHomds.ncoords+(1:3);
HTHomds.phases = 1:HTHomds.nphase;
HTHomds.ntstcol = HTHomds.ntst*HTHomds.ncol;
HTHomds.wp = kron(HTHomds.wpvec',eye(HTHomds.nphase));
HTHomds.pwwt = kron(HTHomds.wt',eye(HTHomds.nphase));
HTHomds.pwi = HTHomds.wi(ones(1,HTHomds.nphase),:);

HTHomds.bialt_M1 = [];
HTHomds.bialt_M2 = [];
HTHomds.bialt_M3 = [];
HTHomds.bialt_M4 = [];
HTHomds.multipliers = nan;
HTHomds.monodromy = [];
HTHomds.multi_r1 = [];
HTHomds.multi_r2 = [];
HTHomds.ups = [];
HTHomds.vps = [];
HTHomds.tsts = 1:HTHomds.ntst;
HTHomds.cols = 1:HTHomds.ncol;