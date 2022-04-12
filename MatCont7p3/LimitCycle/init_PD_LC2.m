function [x,v] = init_PD_LC2(odefile, x, s, ap, ntst, ncol, h)
%
% [x0,v0] = init_PD_LC(odefile, x, s, ntst, ncol)
%
% Initializes a limit cycle continuation of a double period cycle
% from a Period Doubling bifurcation detected during a previous
% limit cycle continuation.
%
global lds cds
% initialize lds
init_lds(odefile,x,s,ap,ntst,ncol);
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


lds.odefile = odefile;
lds.func = func_handles{2};
lds.Jacobian  = func_handles{3};
lds.JacobianP = func_handles{4};
lds.Hessians  = func_handles{5};
lds.HessiansP = func_handles{6};
lds.Der3 = func_handles{7};
lds.Der4 = func_handles{8};
lds.Der5 = func_handles{9};

siz = size(func_handles,2);
if siz > 9
    j=1;
    for i=10:siz
        lds.user{j}= func_handles{i};
        j=j+1;
    end
else
    lds.user=[];
end


lds.nphase = floor((size(x,1)-2)/(s.data.ntst*s.data.ncol+1));
lds.P0 = s.data.parametervalues;
lds.T = 2*s.data.T;
set_ntst_ncol(s.data.ntst,s.data.ncol,s.data.timemesh);

% get x
x = x(:,s.index);
n_par = size(lds.ActiveParams,2);
if n_par==2
     x=[x(lds.coords);lds.P0(lds.ActiveParams)];
else
    x=[x(lds.coords);lds.T;lds.P0(lds.ActiveParams)];
end

% generate double period solution
compute_phi(x);
s.data.phi=lds.PD_phi;
w = s.data.phi/norm(s.data.phi);
w=w(:)';

 if n_par==1
    x = [x(1:(lds.ncoords-lds.nphase))+h*w(1:(lds.ncoords-lds.nphase))';x(lds.coords)-h*w';lds.T;x(end)];
 else
     x = [x(1:(lds.ncoords-lds.nphase))+h*w(1:(lds.ncoords-lds.nphase))';x(lds.coords)-h*w';lds.P0(lds.ActiveParams)];     
 end
 
v = [w(1:(lds.ncoords-lds.nphase))';-w';0;0];
set_ntst_ncol(2*s.data.ntst,s.data.ncol,[lds.msh./2 (lds.msh(2:end)+1)./2])

% generate a new mesh and interpolate
[x,v]=new_mesh(x,v,ntst,ncol);

function init_lds(odefile,x,s,ap,ntst,ncol)
global lds
lds=[];
lds.odefile = odefile;
func_handles = feval(lds.odefile);
lds.func = func_handles{2};
lds.Jacobian  = func_handles{3};
lds.JacobianP = func_handles{4};
lds.Hessians  = func_handles{5};
lds.HessiansP = func_handles{6};
lds.Der3=func_handles{7};
lds.Der4=func_handles{8};
lds.Der5=func_handles{9};
siz = size(func_handles,2);
if siz > 9
    j=1;
    for k=10:siz
        lds.user{j}= func_handles{k};
        j=j+1;
    end
else lds.user=[];
end
lds.nphase = floor((size(x,1)-2)/(s.data.ntst*s.data.ncol+1));
lds.ActiveParams = ap;
lds.P0 = s.data.parametervalues;
set_ntst_ncol(s.data.ntst,s.data.ncol,s.data.timemesh);
lds.T = [];
lds.cols_p1 = 1:(lds.ncol+1);
lds.cols_p1_coords = 1:(lds.ncol+1)*lds.nphase;
lds.ncol_coord = lds.ncol*lds.nphase;
lds.col_coords = 1:lds.ncol*lds.nphase;
lds.pars = lds.ncoords+(1:3);
lds.phases = 1:lds.nphase;
lds.ntstcol = lds.ntst*lds.ncol;
lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwwt = kron(lds.wt',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);

lds.PD_psi = [];
lds.PD_phi = [];
lds.PD_new_phi = [];
lds.PD_new_psi = [];
lds.PD_switch = 0;

lds.BP_psi = [];
lds.BP_phi = [];
lds.BP_psi1 = [];
lds.BP_phi1 = [];
lds.BP_new_phi = [];
lds.BP_new_psi = [];
lds.BP_new_psi1 = [];
lds.BP_new_phi1 = [];
lds.BP_switch = 0;

lds.BPC_switch = 0;
lds.BPC_psi = [];
lds.BPC_phi1 = [];
lds.BPC_phi2 = [];

lds.LPC_phi = [];
lds.LPC_psi = [];
lds.LPC_new_phi = [];
lds.LPC_new_psi = [];
lds.LPC_switch = 0;

lds.NS_psi0 = [];
lds.NS_psi1 = [];
lds.NS_phi0 = [];
lds.NS_phi1 = [];
lds.NS1_new_phi = [];
lds.NS2_new_phi = [];
lds.NS1_new_psi = [];
lds.NS2_new_psi = [];
lds.NS_new_phi = [];
lds.NS_new_psi = [];
lds.NS_switch = 0;
lds.NS1_switch = 0;
lds.NS2_switch = 0;


lds.bialt_M1 = [];
lds.bialt_M2 = [];
lds.bialt_M3 = [];
lds.bialt_M4 = [];
lds.multipliers = nan;
lds.monodromy = [];
lds.multi_r1 = [];
lds.multi_r2 = [];
lds.BranchParam = lds.ActiveParams;
lds.ups = [];
lds.vps = [];
lds.BranchParams=[];
lds.tsts = 1:lds.ntst;
lds.cols = 1:lds.ncol;

function compute_phi(x)
global lds
[xt,p,T] = rearr(x);
ups = reshape(xt,lds.nphase,lds.tps);
p = num2cell(p);

jac = spalloc(lds.ncoords,lds.ncoords,(lds.ncol+4)*lds.nphase);

range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
for j=1:lds.ntst
  xp = ups(:,range0)*lds.wt;
  jac(range1,range2) = bordBVP_PD_jac_f(lds.func,xp,p,T,j);

  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end

% boundary conditions
range = (lds.tps-1)*lds.nphase + (lds.phases);
jac(range,[lds.phases range]) = bordBVP_PD_jac_bc;

%compute borders
jace=[jac rand(lds.ncoords,1);rand(1,lds.ncoords) 0];
 b = []; b(lds.ncoords+1)=1;
  q=jace\b'; q=q(1:end-1);q=q/norm(q);
  p=jace'\b';p=p(1:end-1);p=p/norm(p);
  lds.PD_phi = reshape(q,lds.nphase,lds.tps);
  lds.PD_psi = reshape(p,lds.nphase,lds.tps);
  
function [x,p,T] = rearr(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x), parameters (p) and period (T)
global lds

nap = length(lds.ActiveParams);

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+(1:nap));
x = x0(lds.coords);
T = x0(lds.PeriodIdx);
lds.T = T;

