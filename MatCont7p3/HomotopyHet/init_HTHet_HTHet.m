function [x,v] = init_HTHet_HTHet(odefile, x, v, s, p, ap, up, aup, sp, asp, ntst, ncol,T,eps1,eps1tol)

global HTHetds cds

oldHTHetds = HTHetds;
% initialize HTHetds

[y,v] = init_HTHetds(odefile,x,v,p,ap,up,aup,sp,asp,ntst,ncol,T,eps1,eps1tol,s);

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


HTHetds.odefile = odefile;
HTHetds.func = func_handles{2};
HTHetds.Jacobian  = func_handles{3};
HTHetds.JacobianP = func_handles{4};
HTHetds.Hessians  = func_handles{5};
HTHetds.HessiansP = func_handles{6};
HTHetds.Der3 = func_handles{7};
HTHetds.Der4 = func_handles{8};
HTHetds.Der5 = func_handles{9};


A1 = cjac(HTHetds.func,HTHetds.Jacobian,HTHetds.x0,num2cell(p),HTHetds.ActiveParams);
D1 = eig(A1);
A2 = cjac(HTHetds.func,HTHetds.Jacobian,HTHetds.x1,num2cell(p),HTHetds.ActiveParams);
D2 = eig(A2);

% nneg = dimension of stable subspace
HTHetds.npos = sum(real(D1) > 0);
if (HTHetds.npos == HTHetds.nphase)
    if min(abs(real(D1))) < 1e-10
        HTHetds.npos = HTHetds.npos -1;
    end
end
if (HTHetds.npos == 0)
    if min(abs(real(D1))) < 1e-10
        HTHetds.npos = HTHetds.npos +1;
    end
end

HTHetds.nneg = sum(real(D2) < 0);
if (HTHetds.nneg == HTHetds.nphase)
    if min(abs(real(D2))) < 1e-10
        HTHetds.nneg = HTHetds.nneg -1;
    end
end
if (HTHetds.nneg == 0)
    if min(abs(real(D2))) < 1e-10
        HTHetds.nneg = HTHetds.nneg +1;
    end
end


if HTHetds.index == 1
    if ~((length(HTHetds.ActiveUParams)+length(HTHetds.ActiveSParams))==(HTHetds.nphase-HTHetds.nneg+2))
        message = strcat('The sum of free connection parameters should be equal to ',' ',num2str(HTHetds.nphase-HTHetds.nneg+2));
        error(message);
    end
elseif HTHetds.index == 2
    if ~((length(HTHetds.ActiveParams)+length(HTHetds.ActiveSParams))==(HTHetds.nphase-HTHetds.nneg-HTHetds.npos+2))
        message = strcat('The sum of free stable parameters and free parameters should be equal to ',' ',num2str(HTHetds.nphase-HTHetds.nneg-HTHetds.npos+2));
        error(message);
    end
else
    if ~((length(HTHetds.ActiveParams))==(HTHetds.nphase-HTHetds.nneg-HTHetds.npos+1))
        message = strcat(num2str(HTHetds.nphase-HTHetds.nneg-HTHetds.npos+1),' parameters should be denoted as free');
        error(message);
    end
end

HTHetds.YUsize = (HTHetds.nphase-HTHetds.npos)*HTHetds.npos;
HTHetds.YSsize = (HTHetds.nphase-HTHetds.nneg)*HTHetds.nneg;
% COMPOSE X0
% % ----------
if HTHetds.index == 1
    
    % 1. cycle 
    ups = reshape(y,HTHetds.nphase,HTHetds.tps);
    HTHetds.upold = ups;

    % 2. free UParams
    y = [y; HTHetds.UParams(HTHetds.ActiveUParams)];
    % 3. free SParams
    y = [y; HTHetds.SParams(HTHetds.ActiveSParams)];
    % 4. eps1
    y = [y; HTHetds.eps1];
    
    HTHetds.PeriodIdx = HTHetds.ncoords+2*HTHetds.nphase;
elseif HTHetds.index == 2
    
    % 1. cycle 
    ups = reshape(y,HTHetds.nphase,HTHetds.tps);
    HTHetds.upold = ups;

    % 2. equilibrium coordinates
    y = [y; HTHetds.x0;HTHetds.x1];
    HTHetds.PeriodIdx = length(y);
    % 3. free parameters
    y = [y; HTHetds.P0(ap)];
    % 4. free SParams
    y = [y; HTHetds.SParams(HTHetds.ActiveSParams)];
    % 5. eps1
    y = [y; HTHetds.eps1];
    % 6. YS and YU, initialized to 0
    %eerst YU, dan YS
    for i=1:HTHetds.nphase-HTHetds.npos
        y = [y; zeros(HTHetds.npos,1)];
    end
    for i=1:HTHetds.nphase-HTHetds.nneg
        y = [y; zeros(HTHetds.nneg,1)];
    end

else
    
    % 1. cycle 
    ups = reshape(y,HTHetds.nphase,HTHetds.tps);
    HTHetds.upold = ups;

    % 2. equilibrium coordinates
    y = [y; HTHetds.x0;HTHetds.x1];
    HTHetds.PeriodIdx = length(y);
    % 3. free parameters
    y = [y; HTHetds.P0(ap)];
    % 4. T
    y = [y; HTHetds.T];
    % 5. eps1
    y = [y; HTHetds.eps1];
    % 6. YS and YU, initialized to 0
    for i=1:HTHetds.nphase-HTHetds.npos
        y = [y; zeros(HTHetds.npos,1)];
    end
    for i=1:HTHetds.nphase-HTHetds.nneg
        y = [y; zeros(HTHetds.nneg,1)];
    end
    
end

x = y;
v = [];

% ASSIGN SOME VALUES TO HOMOCLINIC FIELDS
% ---------------------------------------

HTHetds.YS = zeros(HTHetds.nphase-HTHetds.nneg,HTHetds.nneg);
HTHetds.YU = zeros(HTHetds.nphase-HTHetds.npos,HTHetds.npos);

[QS, eigvlS] = computeBaseHTHet(A2,1,HTHetds.nneg);
[QU, eigvlU] = computeBaseHTHet(A1,0,HTHetds.npos);

HTHetds.oldStableQ = QS;
HTHetds.oldUnstableQ = QU;

HTHetds.ups = [];
HTHetds.upold = [];
HTHetds.upoldp = [];


%-----------------------------------------------------------------
function [y,v1] = init_HTHetds(odefile,x,v,p,ap,up,aup,sp,asp,ntst,ncol,T,eps1,eps1tol,s)
global HTHetds

oldHTHetds = HTHetds;
HTHetds = [];

HTHetds.odefile = odefile;
func_handles = feval(HTHetds.odefile);
HTHetds.func = func_handles{2};
HTHetds.Jacobian  = func_handles{3};
HTHetds.JacobianP = func_handles{4};
HTHetds.Hessians  = func_handles{5};
HTHetds.HessiansP = func_handles{6};
HTHetds.Der3=[];
siz = size(func_handles,2);
if siz > 9
    j=1;
    for k=10:siz
        HTHetds.user{j}= func_handles{k};
        j=j+1;
    end
else
    HTHetds.user=[];
end

HTHetds.nphase = oldHTHetds.nphase;
HTHetds.ntst = ntst;
HTHetds.ncol = ncol;
HTHetds.tps = ntst*ncol+1;
HTHetds.ncoords = HTHetds.tps*oldHTHetds.nphase;
HTHetds.coords = 1:HTHetds.ncoords;
HTHetds.msh = s.data.timemesh;
HTHetds.x0 = x(oldHTHetds.nphase*(oldHTHetds.ntst*oldHTHetds.ncol+1)+1:oldHTHetds.nphase*(oldHTHetds.ntst*oldHTHetds.ncol+1)+oldHTHetds.nphase);
HTHetds.x1 = x(oldHTHetds.nphase*(oldHTHetds.ntst*oldHTHetds.ncol+1)+oldHTHetds.nphase+1:oldHTHetds.nphase*(oldHTHetds.ntst*oldHTHetds.ncol+1)+2*oldHTHetds.nphase);
HTHetds.ActiveParams = ap;
HTHetds.ActiveUParams = aup;
HTHetds.ActiveSParams = asp;
HTHetds.P0 = p;
HTHetds.UParams = up;
HTHetds.SParams = sp;
HTHetds.index = oldHTHetds.index;

%%%
y = x(1:(oldHTHetds.ntst*oldHTHetds.ncol+1)*oldHTHetds.nphase);

if isempty(v)
    v1 = [];
else
    v1 = v(1:(oldHTHetds.ntst*oldHTHetds.ncol+1)*oldHTHetds.nphase);
end
[y,v1] = HTHet_new_mesh(y,v1,oldHTHetds.ntst,oldHTHetds.ncol);

HTHetds.T = T;
HTHetds.eps0 = oldHTHetds.eps0;
HTHetds.eps1 = eps1;

HTHetds.TestTolerance = oldHTHetds.TestTolerance;
HTHetds.eps1tol = eps1tol;

HTHetds.cols_p1 = 1:(HTHetds.ncol+1);
HTHetds.cols_p1_coords = 1:(HTHetds.ncol+1)*HTHetds.nphase;
HTHetds.ncol_coord = HTHetds.ncol*HTHetds.nphase;
HTHetds.col_coords = 1:HTHetds.ncol*HTHetds.nphase;
HTHetds.pars = HTHetds.ncoords+2*HTHetds.nphase+(1:length(ap));
HTHetds.phases = 1:HTHetds.nphase;
HTHetds.ntstcol = HTHetds.ntst*HTHetds.ncol;
HTHetds.wp = kron(HTHetds.wpvec',eye(HTHetds.nphase));
HTHetds.pwwt = kron(HTHetds.wt',eye(HTHetds.nphase));
HTHetds.pwi = HTHetds.wi(ones(1,HTHetds.nphase),:);

HTHetds.bialt_M1 = [];
HTHetds.bialt_M2 = [];
HTHetds.bialt_M3 = [];
HTHetds.bialt_M4 = [];
HTHetds.multipliers = nan;
HTHetds.monodromy = [];
HTHetds.multi_r1 = [];
HTHetds.multi_r2 = [];
HTHetds.ups = [];
HTHetds.vps = [];
HTHetds.tsts = 1:HTHetds.ntst;
HTHetds.cols = 1:HTHetds.ncol;
