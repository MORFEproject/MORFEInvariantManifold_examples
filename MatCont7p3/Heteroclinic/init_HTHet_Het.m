function [x,v] = init_HTHet_Het(odefile, x, v, s, p, ap, ntst, ncol,extravec,T,eps0,eps1)

global hetds cds

oldhetds = hetds;
% initialize hetds

[y,v] = init_hetds(odefile,x,v,p,ap,ntst,ncol,extravec,T,eps0,eps1,s);

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


hetds.odefile = odefile;
hetds.func = func_handles{2};
hetds.Jacobian  = func_handles{3};
hetds.JacobianP = func_handles{4};
hetds.Hessians  = func_handles{5};
hetds.HessiansP = func_handles{6};
hetds.Der3 = func_handles{7};
hetds.Der4 = func_handles{8};
hetds.Der5 = func_handles{9};


A1 = cjac(hetds.func,hetds.Jacobian,hetds.x0,num2cell(p),hetds.ActiveParams);
D1 = eig(A1);
A2 = cjac(hetds.func,hetds.Jacobian,hetds.x1,num2cell(p),hetds.ActiveParams);
D2 = eig(A2);

% nneg = dimension of stable subspace
hetds.npos = sum(real(D1) > 0);
if (hetds.npos == hetds.nphase)
    if min(abs(real(D1))) < 1e-10
        hetds.npos = hetds.npos -1;
    end
end
if (hetds.npos == 0)
    if min(abs(real(D1))) < 1e-10
        hetds.npos = hetds.npos +1;
    end
end

hetds.nneg = sum(real(D2) < 0);
if (hetds.nneg == hetds.nphase)
    if min(abs(real(D2))) < 1e-10
        hetds.nneg = hetds.nneg -1;
    end
end
if (hetds.nneg == 0)
    if min(abs(real(D2))) < 1e-10
        hetds.nneg = hetds.nneg +1;
    end
end

if hetds.npos == 0
    error('It is not possible to start from this point (n0 = 0).');
end
if hetds.nneg == 0
    error('It is not possible to start from this point (n1 = 0).');
end

% check input
n_par = size(ap,2);
if n_par ~= (hetds.nphase-hetds.nneg-hetds.npos+2)
    message = strcat('The sum of free parameters should be equal to ',' ',num2str(hetds.nphase-hetds.nneg-hetds.npos+2));
    error(message);
end

% COMPOSE X0
% % ----------
 
% 1. cycle 
ups = reshape(y,hetds.nphase,hetds.tps);
hetds.upold = ups;

% 2. equilibrium coordinates
y = [y; hetds.x0;hetds.x1];
hetds.PeriodIdx = length(y);
% 3. free parameters
y = [y; hetds.P0(ap)];
% 4. extra free parameters
extravec = [hetds.T; hetds.eps0; hetds.eps1];
y = [y; extravec(find(hetds.extravec))];
% 5. YS and YU, initialized to 0
for i=1:hetds.nphase-hetds.npos
    y = [y; zeros(hetds.npos,1)];
end
for i=1:hetds.nphase-hetds.nneg
    y = [y; zeros(hetds.nneg,1)];
end

x = y;
v = [];


% ASSIGN SOME VALUES TO HOMOCLINIC FIELDS
% ---------------------------------------

hetds.YS = zeros(hetds.nphase-hetds.nneg,hetds.nneg);
hetds.YU = zeros(hetds.nphase-hetds.npos,hetds.npos);
hetds.YUsize = (hetds.nphase-hetds.npos)*hetds.npos;
hetds.YSsize = (hetds.nphase-hetds.nneg)*hetds.nneg;

% Third parameter = unstable_flag, 
% 1 if we want the unstable space, 0 if we want the stable one
[QS, eigvlS] = computeBaseHet(A2,1,hetds.nneg);
[QU, eigvlU] = computeBaseHet(A1,0,hetds.npos);

hetds.oldStableQ = QS;
hetds.oldUnstableQ = QU;

hetds.ups = [];
hetds.upold = [];
hetds.upoldp = [];



%-----------------------------------------------------------------
function [y,v1] = init_hetds(odefile,x,v,p,ap,ntst,ncol,extravec,T,eps0,eps1,s)
global hetds HTHetds
oldhetds = hetds;
hetds = [];

hetds.odefile = odefile;
func_handles = feval(hetds.odefile);
hetds.func = func_handles{2};
hetds.Jacobian  = func_handles{3};
hetds.JacobianP = func_handles{4};
hetds.Hessians  = func_handles{5};
hetds.HessiansP = func_handles{6};
hetds.Der3=[];
siz = size(func_handles,2);
if siz > 9
    j=1;
    for k=10:siz
        hetds.user{j}= func_handles{k};
        j=j+1;
    end
else
    hetds.user=[];
end

hetds.nphase = HTHetds.nphase;
hetds.ntst = HTHetds.ntst;
hetds.ncol = HTHetds.ncol;
hetds.tps = HTHetds.ntst*HTHetds.ncol+1;
hetds.ncoords = HTHetds.tps*HTHetds.nphase;
hetds.coords = 1:HTHetds.ncoords;
hetds.msh = HTHetds.msh;

hetds.x0 = x(end-2*hetds.nphase+1:end-hetds.nphase);
hetds.x1 = x(end-hetds.nphase+1:end);
hetds.ActiveParams = ap;
hetds.P0 = p;
hetds.extravec = extravec;


y = x(1:hetds.ncoords);
if isempty(v)
    v1 = [];
else
    v1 = v(1:hetds.ncoords);
end

[y,v1] = Het_new_mesh(y,v1,ntst,ncol);

hetds.T = T;
hetds.eps0 = eps0;
hetds.eps1 = eps1;

hetds.cols_p1 = 1:(hetds.ncol+1);
hetds.cols_p1_coords = 1:(hetds.ncol+1)*hetds.nphase;
hetds.ncol_coord = hetds.ncol*hetds.nphase;
hetds.col_coords = 1:hetds.ncol*hetds.nphase;
hetds.pars = hetds.ncoords+2*hetds.nphase+(1:length(ap));
hetds.phases = 1:hetds.nphase;
hetds.ntstcol = hetds.ntst*hetds.ncol;
hetds.wp = kron(hetds.wpvec',eye(hetds.nphase));
hetds.pwwt = kron(hetds.wt',eye(hetds.nphase));
hetds.pwi = hetds.wi(ones(1,hetds.nphase),:);

hetds.bialt_M1 = [];
hetds.bialt_M2 = [];
hetds.bialt_M3 = [];
hetds.bialt_M4 = [];
hetds.multipliers = nan;
hetds.monodromy = [];
hetds.multi_r1 = [];
hetds.multi_r2 = [];
hetds.ups = [];
hetds.vps = [];
hetds.tsts = 1:hetds.ntst;
hetds.cols = 1:hetds.ncol;
