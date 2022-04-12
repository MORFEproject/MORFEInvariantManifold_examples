function [x,v] = init_HTHom_Hom(odefile, x, v, s, p, ap, ntst, ncol,extravec,T,eps0,eps1)
%maar 1 rij van x en v ingeven!

global homds cds

% check input
n_par = size(ap,2);
if n_par ~= 2
    error('2 free system parameters are needed');
end

oldhomds = homds;

% initialize homds
[x1,v] = init_homds(odefile,x,v,p,ap,ntst,ncol,extravec,T,eps0,eps1,s);

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


homds.odefile = odefile;
homds.func = func_handles{2};
homds.Jacobian  = func_handles{3};
homds.JacobianP = func_handles{4};
homds.Hessians  = func_handles{5};
homds.HessiansP = func_handles{6};
homds.Der3 = func_handles{7};
homds.Der4 = func_handles{8};
homds.Der5 = func_handles{9};



A = cjac(homds.func,homds.Jacobian,homds.x0,num2cell(p),homds.ActiveParams);
D = eig(A);
% nneg = dimension of stable subspace
homds.nneg = sum(real(D) < 0);
if (homds.nneg == homds.nphase)
    if min(abs(real(D))) < 1e-2
        homds.nneg = homds.nneg -1;
    end
end
if (homds.nneg == 0)
    if min(abs(real(D))) < 1e-2
        homds.nneg = homds.nneg +1;
    end
end
homds.npos = homds.nphase-homds.nneg;
homds.Ysize = homds.nneg*homds.npos;
if homds.nneg == 0
    error('It is not possible to start from this point (nS = 0)');
end
if homds.npos == 0
    error('It is not possible to start from this point (nU = 0)');
end

% COMPOSE X0
% % ----------

% 1. cycle 
ups = reshape(x1,homds.nphase,homds.tps);
homds.upold = ups;

% 2. equilibrium coordinates
x1 = [x1; homds.x0];
homds.PeriodIdx = length(x1);
% 3. (two) free parameters
x1 = [x1; homds.P0(ap)];
% 4. extra free parameters
extravec = [homds.T; homds.eps0; homds.eps1];
x1 = [x1; extravec(find(homds.extravec))];
% 5. YS and YU, initialized to 0
for i=1:homds.nneg
    x1 = [x1; zeros(homds.npos,1)];
end
for i=1:homds.npos
    x1 = [x1; zeros(homds.nneg,1)];
end

x = x1;
v = [];


% ASSIGN SOME VALUES TO HOMOCLINIC FIELDS
% ---------------------------------------

homds.YS = zeros(homds.npos,homds.nneg);
homds.YU = zeros(homds.nneg,homds.npos);

% Third parameter = unstable_flag, 
% 1 if we want the unstable space, 0 if we want the stable one
[QU, se] = computeBase(A,0,homds.npos);
[QS, se] = computeBase(A,1,homds.nneg);

homds.oldStableQ = QS;
homds.oldUnstableQ = QU;
homds.ups = [];
homds.upold = [];
homds.upoldp = [];

%-----------------------------------------------------------------
function [x1,v] = init_homds(odefile,x,v,p,ap,ntst, ncol,extravec,T,eps0,eps1,s)
global homds HTHomds
oldhomds = homds;
homds = [];

homds.odefile = odefile;
func_handles = feval(homds.odefile);
homds.func = func_handles{2};
homds.Jacobian  = func_handles{3};
homds.JacobianP = func_handles{4};
homds.Hessians  = func_handles{5};
homds.HessiansP = func_handles{6};
homds.Der3=[];
siz = size(func_handles,2);
if siz > 9
    j=1;
    for k=10:siz
        homds.user{j}= func_handles{k};
        j=j+1;
    end
else homds.user=[];end

homds.nphase = HTHomds.nphase;
homds.x0 = x(HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.nphase);
homds.ActiveParams = ap;
homds.P0 = p;
homds.extravec = extravec;

%%%
homds.ntst = HTHomds.ntst;
homds.ncol = HTHomds.ncol;
homds.tps = HTHomds.ntst*HTHomds.ncol+1;
homds.ncoords = HTHomds.tps*HTHomds.nphase;
homds.coords = 1:HTHomds.ncoords;
homds.msh = HTHomds.msh;

x1 = x(1:homds.ncoords);
[x1,v]=Hom_new_mesh(x1,v,ntst,ncol);


homds.T = T;
homds.eps0 = eps0;
homds.eps1 = eps1;
homds.cols_p1 = 1:(homds.ncol+1);
homds.cols_p1_coords = 1:(homds.ncol+1)*homds.nphase;
homds.ncol_coord = homds.ncol*homds.nphase;
homds.col_coords = 1:homds.ncol*homds.nphase;
homds.pars = homds.ncoords+(1:3);
homds.phases = 1:homds.nphase;
homds.ntstcol = homds.ntst*homds.ncol;
homds.wp = kron(homds.wpvec',eye(homds.nphase));
homds.pwwt = kron(homds.wt',eye(homds.nphase));
homds.pwi = homds.wi(ones(1,homds.nphase),:);

homds.bialt_M1 = [];
homds.bialt_M2 = [];
homds.bialt_M3 = [];
homds.bialt_M4 = [];
homds.multipliers = nan;
homds.monodromy = [];
homds.multi_r1 = [];
homds.multi_r2 = [];
homds.ups = [];
homds.vps = [];
homds.tsts = 1:homds.ntst;
homds.cols = 1:homds.ncol;

homds.HTPstep = 0;