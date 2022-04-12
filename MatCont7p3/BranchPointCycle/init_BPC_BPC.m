function [x,v] = init_BPC_BPC(odefile, x, s, ap, ntst, ncol,bp)
%
% [x0,v0] = init_BPC_BPC(odefile, x, s, ap, ntst, ncol, bp)
%
global lds cds

% check input
if size(ap,2)~=3
  error('Two active parameters are needed for a Branch Point bifurcation curve continuation');
end
if size(bp,2)~=1
    errordlg('One branch parameter is needed for a Branchpoint bifurcation curve continuation');
end

% initialize lds
init_lds(odefile,x,s,ap,ntst,ncol,bp);

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

% get x
x = x(:,s.index);

% append extra parameter
x(lds.ncoords+(2:4),:) = lds.P0(ap);

% generate a new mesh and interpolate
[x]=new_mesh(x,[],ntst,ncol);

% let continuer calculate tangent vector
v = [];

[xt,p,T] = rearr(x);
ups = reshape(xt,lds.nphase,lds.tps);
p = num2cell(p);
pars1 = lds.ncoords+1;
pars2 = lds.ncoords+2;
jac = spalloc(lds.ncoords+1,lds.ncoords+2,(lds.ncol+4)*lds.nphase);

% function
range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
for j=lds.tsts
  xp = ups(:,range0)*lds.wt;
  jac(range1,[range2 pars1 pars2]) = bordBVP_BPC_f(lds.func,xp,p,T,j);
  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end
% boundary conditions
range  = (lds.tps-1)*lds.nphase+ (lds.phases);
range1 = lds.ncoords-lds.nphase+lds.phases;
jac(range,[lds.phases range1]) = bordBVP_LPC_bc1;
% integral constraint
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1;
range2 = lds.cols_p1_coords;
for j=lds.tsts
  p = lds.dt(j)*(ups(:,range1).*lds.pwi);
  ic(range2) = ic(range2)+p(lds.cols_p1_coords);
  range1 = range1 + lds.ncol;
  range2 = range2 + lds.ncol_coord;
end

jac(range(end)+1,1:lds.ncoords)= ic;

%compute borders
[Q,R,E] = qr(full(jac));
R(end,end)=0;R(end,end-1)=0;
p = E*[R(1:end-1,1:end-2)\-R(1:end-1,end-1:end);eye(2)];
p = p'/norm(p);
lds.BPC_new_phi = [];
lds.BPC_new_psi = [];
lds.BPC_phi = p;
lds.BPC_phi1=p(1,:);
lds.BPC_phi2=p(2,:);
p = Q(:,end);
lds.BPC_psi = p';


% ---------------------------------------------------------------
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

%-----------------------------------------------------------------
function init_lds(odefile,x,s,ap,ntst,ncol,bp)
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
lds.BranchParam = bp;
lds.ups = [];
lds.vps = [];
lds.BranchParams=[];
lds.tsts = 1:lds.ntst;
lds.cols = 1:lds.ncol;
