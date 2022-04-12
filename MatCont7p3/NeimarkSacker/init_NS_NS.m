function [x0,v0] = init_NS_NS(odefile, x, s, ap, ntst, ncol)
%
% [x0,v0] = init_NS_NS(odefile, x, s, ap, ntst, ncol)
%
global lds cds
% check input
if size(ap)~=[1 2]
  error('Two active parameters are needed for a torus bifurcation curve continuation');
end
% initialize lds
% This will (among other things) correctly set the parameters (fixed and free) and T in lds.
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

x1 = x(1:lds.ncoords,s.index);
p2 = lds.P0;
p1 = num2cell(p2);
T = lds.T;
% append period and extra parameters
x0 = x1;
x0(lds.ncoords+1) = lds.T;
x0(lds.ncoords+(2:3)) = lds.P0(ap);

% generate a new mesh and interpolate
x0 = new_mesh(x0,[],ntst,ncol);
x1 = x0(1:lds.ncoords);
v0 = [];
%find pair of complex eigenvalues with modulus 1
if (~isfield(s.data,'multipliers') || isempty(s.data.multipliers))
    error('Multipliers not provided in s.data')
end
d = s.data.multipliers;
smallest_mod = Inf;
idx1=[];
for j=1:lds.nphase
    val=find(d==conj(d(j)));
    if ~isempty(val) && val~=j
        idx=abs(abs(d(val))-1);
        if idx < smallest_mod
            idx1 = j;
            idx2 = val;
            smallest_mod = idx;
        end
    end
end
smallest_sum = Inf;
if isempty(idx1)
    for j=1:lds.nphase-1
        [val,idx] = min(abs(d(j+1:lds.nphase)*d(j)-1));
        if val < smallest_sum
            idx1 = j;
            idx2 = j+idx;
            smallest_sum = val;
        end
    end
end
k = (d(idx2)+d(idx1))/2;

% append k
x0(lds.ncoords+4) = k;

jac = spalloc(2*lds.ncoords-lds.nphase,2*lds.ncoords-lds.nphase,(2*lds.tps-1)*lds.nphase*3+(2*lds.ncol+1)*lds.nphase*2*lds.ntst);
ups = reshape(x1,lds.nphase,lds.tps);
% function

range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
for ns = 1:2
    range0 = lds.cols_p1;
    for j = lds.tsts
        xp = ups(:,range0)*lds.wt;
        jac(range1,range2) = bordBVP_NS_f(lds.func,xp,p1,T,j);
        range0 = range0 + lds.ncol;
        range1 = range1 + lds.ncol_coord;
        range2 = range2 + lds.ncol_coord;
    end
end
% boundary conditions
range  = 2*(lds.tps-1)*lds.nphase+ (lds.phases);
range2  = 2*(lds.ncoords-lds.nphase)+lds.phases;
range1 = lds.ncoords-lds.nphase+lds.phases;
jac(range,[lds.phases range1 range2]) = bordBVP_NS_bc1(k);

dimdoub=2*lds.ncoords-lds.nphase;
jace=[jac rand(dimdoub,2);rand(2,dimdoub) zeros(2,2)];
b=zeros(dimdoub+2,2);b(dimdoub+1,1)=1;b(dimdoub+2,2)=1;
q=jace\b;q=q(1:dimdoub,:);q=orth(q);
lds.NS_phi0 = q(:,1); 
lds.NS_phi1 = q(:,2);
p=jace'\b;p=p(1:dimdoub,:);p=orth(p);
lds.NS_psi0 = p(:,1); 
lds.NS_psi1 = p(:,2);

for i=1:lds.tps
    lds.upoldp(:,i) = T*feval(lds.func, 0, ups(:,i), p1{:});
end

%compute indices

A = NS_BVP_jac('BVP_LC_jac_f','BVP_LC_jac_bc','BVP_LC_jac_ic',x1,p2,T,3,2);
[Q,R] = qr(A');
Bord  = [jac lds.NS_psi0 lds.NS_psi1; [lds.NS_phi0';lds.NS_phi1'] zeros(2)];
bunit = zeros((2*lds.tps-1)*lds.nphase+2,2);
bunit(end-1:end,:) = eye(2);
sn = Bord\bunit;
st = Bord'\bunit;
v = sn(1:end-2,:)';
w = st(1:end-2,:)';
w1 = w(:,end-lds.nphase+1:end);
% calculate g'
ups = reshape(x1,lds.nphase,lds.tps);

range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
cv=[];
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
gx = zeros(4,lds.ncoords+3);gk=[];
for ns = 1:2
      range3 = lds.cols_p1_coords;
      v1     = cv  ;
      range0 = lds.cols_p1;
      for tstpt = lds.tsts
          xp  = ups(:,range0)*lds.wt;
          cv = v(:,range2)';
          cw = w(:,range1);
          range = lds.phases;
          for c = lds.cols
              xt    = xp(:,c);
              sysj  = cjac(lds.func,lds.Jacobian,xt,p1,lds.ActiveParams);
              sysh  = chess(lds.func,lds.Jacobian,lds.Hessians,xt,p1,lds.ActiveParams);
              syshp = chessp(lds.func,lds.Jacobian,lds.HessiansP,xt,p1,lds.ActiveParams);
              wtk   = lds.wt(kr1,c(ones(1,lds.nphase)))';
              for d = lds.phases
                  sh1(:,d) = (wtk.*sysh(:,kr2,d))*cv(:,1);
                  sh2(:,d) = (wtk.*sysh(:,kr2,d))*cv(:,2);
              end      
              t11 = T* wtk.*sh1(:,kr2);
              t21 = T* wtk.*sh2(:,kr2);
              t12 = (wtk.*sysj(:,kr2))*cv(:,1);
              t22 = (wtk.*sysj(:,kr2))*cv(:,2);
              t13 = T* wtk.*syshp(:,kr2,1)* cv(:,1);
              t23 = T* wtk.*syshp(:,kr2,1)* cv(:,2);
              t14 = T* wtk.*syshp(:,kr2,2)* cv(:,1);
              t24 = T* wtk.*syshp(:,kr2,2)* cv(:,2);
              syshess1(range,:) = [t11 t12 t13 t14];      
              syshess2(range,:) = [t21 t22 t23 t24];      
              range = range + lds.nphase;    
          end
          gx(1,[range3 lds.ncoords+(1:3)]) = gx(1,[range3 lds.ncoords+(1:3)]) + cw(1,:)*syshess1;
          gx(2,[range3 lds.ncoords+(1:3)]) = gx(2,[range3 lds.ncoords+(1:3)]) + cw(1,:)*syshess2;
          gx(3,[range3 lds.ncoords+(1:3)]) = gx(3,[range3 lds.ncoords+(1:3)]) + cw(2,:)*syshess1;
          gx(4,[range3 lds.ncoords+(1:3)]) = gx(4,[range3 lds.ncoords+(1:3)]) + cw(2,:)*syshess2;
          range0 = range0 + lds.ncol;
          range1 = range1 + lds.ncol_coord;
          range2 = range2 + lds.ncol_coord;
          range3 = range3 + lds.ncol_coord;
      end
end  
gk(1,1) = 2*w1(1,:)*v1(end-lds.nphase+1:end,1);
gk(2,1) = 2*w1(1,:)*v1(end-lds.nphase+1:end,2);
gk(3,1) = 2*w1(2,:)*v1(end-lds.nphase+1:end,1);
gk(4,1) = 2*w1(2,:)*v1(end-lds.nphase+1:end,2);
B = [A ; gx gk]*Q;
Jres = B(2+lds.ncoords:end,2+lds.ncoords:end)';
[Q,R,E] = qr(full(Jres));
index = [1 1;1 2;2 1;2 2];
[I,J] = find(E(:,1:2));
lds.index1 = index(I(J(1)),:);
lds.index2 = index(I(J(2)),:);

%-----------------------------------------------------------------
function init_lds(odefile,x,s,ap,ntst,ncol)
global lds
oldT = lds.T;
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
if isfield(s.data,'T')
    lds.T = s.data.T;
else
    lds.T = oldT;
end
lds.cols_p1 = 1:(lds.ncol+1);
lds.cols_p1_coords = 1:(lds.ncol+1)*lds.nphase;
lds.ncol_coord = lds.ncol*lds.nphase;
lds.col_coords = 1:lds.ncol*lds.nphase;
lds.pars = lds.ncoords+(1:4);
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
