function [x0,v0] = init_H_LC(odefile, x, p, ap, h, ntst, ncol)

%
% [x0,v0] = init_H_LC(odefile,x,p,ap,h,ntst,ncol)
%
% Initializes a limit cycle continuation from a Hopf point
% 
%
global cds lds eds hds 

% check input
n_par = size(ap,2);
if n_par~=1&& n_par~= 2
    error('One active parameter and the period or 2 active parameters are needed for limit cycle continuation');
end
lds = [];
% initialize lds
if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end
cds.curve = @equilibrium;
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_jacobian = curvehandles{4};
cds.curve_hessians = curvehandles{5};
func_handles = feval(odefile);
eds.func = func_handles{2};
eds.Jacobian  = func_handles{3};
eds.JacobianP = func_handles{4};
if ~isfield(eds,'ActiveParams')
    init_EP_EP(odefile,x,p,ap);
end

init_lds(odefile,x,p,ap,ntst,ncol);

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
cds.symjac = 1;
cds.symhess = 0;


lds.P0 = p;
cds.oldJac = [];
cds.oldJacX = [];
xp = [x;p(ap)];

% xp = [x;p(ap)];
% check parameters

%% starter...

% compute/extract bifurcation data
cds.ndim = length(xp);
p = num2cell(p);
if ~isempty(hds)
    A = cjac(hds.func,hds.Jacobian,x,p,hds.ActiveParams);
else
    A = cjac(cds.curve_func,cds.curve_jacobian,xp,[]);
end

% nphase = size(x);
nphase = lds.nphase;
A = A(1:nphase,1:nphase);
% calculate eigenvalues and eigenvectors
[V,D] = eig(A);
% find pair of complex eigenvalues
d = diag(D);
smallest_sum = Inf;
for j=1:nphase-1
  [val,idx] = min(abs(d(j+1:nphase)+d(j)));
  if val < smallest_sum
    idx1 = j;
    idx2 = j+idx;
    smallest_sum = val;
  end
end
% real part? Oh dear, a neutral saddle!
if imag(d(idx1)) == 0 & imag(d(idx2)) == 0
    x0=[];v0=[];
  debug('Neutral saddle\n');
  return;
end
% get imaginary part and corresponding eigenvector
omega = abs(imag(d(idx2)));
Q = V(:,idx1);

d = real(Q)'*real(Q);
s = imag(Q)'*imag(Q);
r = real(Q)'*imag(Q);
Q = Q*exp(i*atan2(2*r,s-d)/2);
Q = Q/norm(real(Q));

% initial amplitude h
% calculate initial cycle and its tangent vector
t = kron(exp(2*pi*i*lds.finemsh),Q);
lds.upoldp = -imag(t);
v0 = [real(t(:));0;0];
if n_par==1
    x0 = [repmat(x,lds.tps,1);2*pi/omega;p{lds.ActiveParams}]+h*v0;
else    
    x0 = [repmat(x,lds.tps,1);p{lds.ActiveParams(1)};p{lds.ActiveParams(2)}]+h*v0;
end
lds.T = 2*pi/omega;
v0 = v0/norm(v0);


%-----------------------------------------------------------------
function init_lds(odefile,x,p,ap,ntst,ncol)
global lds cds
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
symord = 0; 
if ~isempty(func_handles{9}),   symord = 5; 
elseif ~isempty(func_handles{8}),   symord = 4; 
elseif ~isempty(func_handles{7}),   symord = 3; 
elseif ~isempty(func_handles{5}),   symord = 2; 
elseif ~isempty(func_handles{3}),   symord = 1; 
end
cds.options = contset(cds.options, 'SymDerivative', symord);
siz = size(func_handles,2);
if siz > 9
    j=1;
    for k=10:siz
        lds.user{j}= func_handles{k};
        j=j+1;
    end
else lds.user=[];
end
lds.nphase = size(x,1);
lds.ActiveParams = ap;
lds.P0 = p;
set_ntst_ncol(ntst,ncol,(0:ntst)/ntst);
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

