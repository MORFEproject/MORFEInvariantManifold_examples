function [x,v] = init_LC_HSN(odefile, x, s, p, ap, ntst, ncol,extravec,T,eps0,eps1)

global homds lds cds

% check input
n_par = size(ap,2);
nextra = length(find(extravec));
if (n_par + nextra) ~= 4
    error('4 active and free homoclinic parameters are needed');
end

if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end
cds.curve = @homoclinicsaddlenode;
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_defaultprocessor = curvehandles{2};
cds.curve_jacobian = curvehandles{4};
cds.curve_hessians = curvehandles{5};
cds.curve_testf = curvehandles{6};
cds.curve_adapt = curvehandles{13};
cds.symjac = 1;
cds.symhess = 0;

if size(x,2) > 1
    x = x(:,s.index);
end
homds = [];

% initialize homds
init_homds(odefile,x,p,ap,lds.ntst,lds.ncol,extravec,T,eps0,eps1);

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
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);


homds.odefile = odefile;
homds.func = func_handles{2};
homds.Jacobian  = func_handles{3};
homds.JacobianP = func_handles{4};
homds.Hessians  = func_handles{5};
homds.HessiansP = func_handles{6};
homds.Der3 = func_handles{7};
homds.Der4 = func_handles{8};
homds.Der5 = func_handles{9};

x = x(1:homds.nphase*homds.tps);

cds.oldJac = [];
cds.oldJacX = [];
xp = [x;p(ap)];
cds.ndim = length(xp);



% STEP 2: Finding and storing the equilibrium
% ===========================================
    
% First find place on the old cycle where dx, the change in coordinates,
% has minimal norm
ups = reshape(x,homds.nphase,homds.tps);
homds.upold = ups;
tmpres(length(ups(1,:))) = 0;
pcell = num2cell(p);
for i=1:length(ups(1,:))
    tmpres(i) = norm(feval(homds.func,0,ups(:,i),pcell{:}));
end
[minval,ind] = min(tmpres);
x0 = ups(:,ind);
    

% STEP 3: Rotating cycle
% ======================
% Build new cycle, cycleold, cycleoldder, mesh and finemesh
    
% The interval containing the equilibrium-approximation will be left out
homds.ntst = homds.ntst-1;
homds.tsts = 1:homds.ntst;
homds.tps = homds.ntst*homds.ncol+1;
homds.ncoords = homds.tps*homds.nphase;
homds.coords = 1:homds.ncoords;
    
meshStartInd = ceil(ind / homds.ncol) + 1;
meshEndInd = ceil(ind / homds.ncol);
meshShiftPoint = length(homds.msh);

cycleStartInd = (meshStartInd-1) * homds.ncol + 1;
cycleEndInd = (meshEndInd-1) * homds.ncol + 1;
cycleShiftPoint = length(homds.finemsh);
    
j = 1;
for i=meshStartInd:meshShiftPoint
    hommesh(j) = homds.msh(i) - homds.msh(meshStartInd);
    j = j+1;
end
if (meshShiftPoint-meshStartInd)
    for i=2:meshEndInd
        hommesh(j) = homds.msh(i) + hommesh(meshShiftPoint - meshStartInd+1);
        j = j+1;
    end
else        
    for i=2:meshEndInd
        hommesh(j) = homds.msh(i);
        j = j+1;
    end
end

newcycleold = [];
newcycleoldder = [];
j = 1;
for i=cycleStartInd:cycleShiftPoint
    newcycle(:,j) = ups(:,i);
    newcycleold(:,j) = homds.upold(:,i);
    if ~isempty(homds.upoldp)
        newcycleoldder(:,j) = homds.upoldp(:,i);
    end
    newfinemesh(j) = homds.finemsh(i) - homds.finemsh(cycleStartInd);
    j = j+1;
end
if (cycleShiftPoint - cycleStartInd)
    for i=2:cycleEndInd
        newcycle(:,j) = ups(:,i);
        newcycleold(:,j) = homds.upold(:,i);
        if ~isempty(homds.upoldp)
            newcycleoldder(:,j) = homds.upoldp(:,i);
        end
        newfinemesh(j) = homds.finemsh(i) + newfinemesh(cycleShiftPoint - cycleStartInd+1);
        j = j+1;
    end
else
    for i=2:cycleEndInd
        newcycle(:,j) = ups(:,i);
        newcycleold(:,j) = homds.upold(:,i);
        if ~isempty(homds.upoldp)
            newcycleoldder(:,j) = homds.upoldp(:,i);
        end
        newfinemesh(j) = homds.finemsh(i);
        j = j+1;
    end
end

ups = newcycle;
homds.upold = newcycleold;
homds.upoldp = newcycleoldder;
homds.msh = hommesh;
homds.dt = homds.msh(2:end) - homds.msh(1:end-1);
homds.finemsh = newfinemesh;
homds.v = [];
    
if (x0 == ups(:,end)) | (x0 == ups(:,1))
    x0 = x0 + 1e-4;
end

homds.x0 = x0;
    
Hom_calc_weights;

A = cjac(homds.func,homds.Jacobian,x0,pcell,homds.ActiveParams);

% nneg = dimension of stable subspace
% one eigenvalue should be (approximately) zero
D=eig(A);

[y,i] = min(abs(D));
homds.nneg = sum(real(D) < 0);
homds.npos = sum(real(D) > 0);

if imag(D(i))
   error('Bad startpoint.');
end
   if D(i) < 0
   homds.nneg = homds.nneg-1;
elseif (homds.nneg+homds.npos == homds.nphase)
   homds.npos = homds.npos-1;
end 


homds.Ysize = (homds.nneg+1)*homds.npos + homds.nneg*(homds.npos+1);

scale = 1/homds.msh(end);
homds.msh = homds.msh * scale;
homds.dt = homds.dt * scale;
homds.finemsh = homds.finemsh * scale;
homds.T = homds.T/scale;

% COMPOSE X0
% ----------

% 1. cycle 
x1 = reshape(ups,size(ups,1)*size(ups,2),1);
v = []; 

[x1,v]=Hom_new_mesh(x1,v,ntst,ncol);

% 2. equilibrium coordinates
x1 = [x1; x0];
% 3. (two) free parameters
homds.PeriodIdx = length(x1);
x1 = [x1; homds.P0(homds.ActiveParams)];
% 4. extra free parameters
extravec = [homds.T; homds.eps0; homds.eps1];
x1 = [x1; extravec(find(homds.extravec))];
% 5. YS and YU, initialized to 0
for i=1:homds.nneg
    x1 = [x1; zeros(homds.npos+1,1)];
end

for i=1:homds.npos
    x1 = [x1; zeros(homds.nneg+1,1)];
end

% ASSIGN SOME VALUES TO HOMOCLINIC FIELDS
% ---------------------------------------

homds.YS = zeros(homds.npos+1,homds.nneg);
homds.YU = zeros(homds.nneg+1,homds.npos);

[Y,i] = min(abs(D));
RED = A - D(i) * eye(homds.nphase);
homds.vvector = null(RED);
homds.wvector = null(RED');

% Third parameter = unstable_flag, 
% 1 if we want the unstable space, 0 if we want the stable one
[QU, eigvlU] = computeBaseHSN(A,0,homds.npos);
[QS, eigvlS] = computeBaseHSN(A,1,homds.nneg);

homds.oldStableQ = QS;
homds.oldUnstableQ = QU;
homds.ups = [];

x = x1;
v = [];



%-----------------------------------------------------------------
function init_homds(odefile,x,p,ap,ntst,ncol,extravec,T,eps0,eps1)
global homds lds
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
else homds.user=[];
end
homds.nphase = lds.nphase;
homds.ActiveParams = ap;
homds.P0 = p;
Hom_set_ntst_ncol(ntst,ncol,(0:ntst)/ntst);
homds.dt = lds.dt;
homds.msh = lds.msh;
homds.finemsh = lds.finemsh;
homds.extravec = extravec;
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

% Border needed for the constraint that the equilibrium is a saddle-node
homds.wvector = [];
homds.vvector = [];
