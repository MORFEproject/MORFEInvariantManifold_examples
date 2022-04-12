function result=nf_FF(x)
%
% calculates FF normal form coefficient
%
global lds

lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);

[xt,p,T] = rearr(x);

pt = num2cell(p);
pars1 = lds.ncoords+1;
jac = spalloc(lds.ncoords,lds.ncoords,lds.ncol*(lds.ncol+1)*lds.nphase^2*lds.ntst);

ups = reshape(xt,lds.nphase,lds.tps);
%compute v1M
% function
range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
for j=lds.tsts
  xp = ups(:,range0)*lds.wt;
  jac(range1,[range2 pars1]) = bordBVP_LPC_f(lds.func,xp,pt,T,j);
  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end
% boundary conditions
range  = (lds.tps-1)*lds.nphase+ (lds.phases);
range1 = lds.ncoords-lds.nphase+lds.phases;
jac(range,[lds.phases range1]) = bordBVP_LPC_bc1;

f1 = -jac(1:lds.ncoords-lds.nphase,lds.ncoords+1);
f1ps = reshape(f1,lds.nphase,lds.tps-1);
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
wt = lds.wt';
for j=lds.tsts
    f1psW1(:,range2) = (f1ps(:,range2).*(repmat(gl_weight(lds.ncol)',lds.nphase,1)/2))*lds.dt(j); %compose f1W1   
    p = f1psW1(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end

jace = jac(lds.coords,lds.coords);
jace = [jace randn(lds.ncoords,1);randn(1,lds.ncoords) 0];
b = []; b(lds.ncoords+1)=1; b=b';
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

jace = [jace(lds.coords,lds.coords) p;ic 0];
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*f1; rl(lds.ncoords+1) = 0;
v1 = jace\rl';
v1 = v1(lds.coords);

jace(lds.ncoords+1,lds.coords) = q';
phi = jace'\b;
phi = phi(1:lds.ncoords-lds.nphase);

%rescale phi
% function
range1 = lds.cols_p1;
range2 = lds.cols;
v1ps = reshape(v1,lds.nphase,lds.tps);
for j=lds.tsts
    % value of polynomial on each collocation point
    v1Cps(:,range2) = v1ps(:,range1)*lds.wt;
    range1 = range1+lds.ncol;
    range2 = range2+lds.ncol;
end
v1C = reshape(v1Cps,lds.ncoords-lds.nphase,1);

ic = phi'*v1C;
phi = phi/ic;

% boundary conditions
Jjac = jace;
range  = (lds.tps-1)*lds.nphase+ (lds.phases);
range1 = lds.ncoords-lds.nphase+lds.phases;
Jjac(range,[lds.phases range1]) = bordBVP_R2_bc1;
Jjac(end,:) = []; Jjac(:,end) = [];

Jjac = [Jjac randn(lds.ncoords,1);randn(1,lds.ncoords) 0];
q=Jjac\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=Jjac'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

Jjac = [Jjac(lds.coords,lds.coords) p;q' 0];
v2 = Jjac\b;
v2 = v2(lds.coords);
v2ps = reshape(v2,lds.nphase,lds.tps);
ficd = dot(v2ps,v2ps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v2 = v2/sqrt(ic);

w2 = Jjac'\b;
w2 = w2(1:lds.ncoords-lds.nphase);

phips = reshape(phi,lds.nphase,lds.tps-1);
v2ps = reshape(v2,lds.nphase,lds.tps);
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
range3 = lds.cols_p1;
wt = lds.wt';
for j=lds.tsts
    v2Cps(:,range2) = v2ps(:,range3)*lds.wt;
    phiCps(:,range2) = (phips(:,range2)./(repmat(gl_weight(lds.ncol)',lds.nphase,1)/2))/lds.dt(j);        
    p = phips(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
    range3 = range3 + lds.ncol;
end
phiC = reshape(phiCps,lds.ncoords-lds.nphase,1);
v2C = reshape(v2Cps,lds.ncoords-lds.nphase,1);

%rescale w2
ic1 = w2'*v2C;
w2 = w2/ic1;

jacee = jace;
jacee(1:lds.ncoords-lds.nphase,lds.ncoords+1) = v1C;
jacee(lds.ncoords-lds.nphase+1:lds.ncoords,lds.ncoords+1) = zeros(lds.nphase,1);
rl = [];
rl(1:lds.ncoords) = T*ic; 
rl(lds.ncoords+1) = 0;
w1 = jacee'\rl';
w1 = w1(1:lds.ncoords-lds.nphase);

% function
range1 = lds.cols_p1;
range2 = lds.phases;
range3 = lds.col_coords;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
w1ps = reshape(w1,lds.nphase,lds.tps-1);
w2ps = reshape(w2,lds.nphase,lds.tps-1);
ic2 = zeros(1,lds.ncoords);
ic3 = zeros(1,lds.ncoords);
wt = lds.wt';
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    v1_1 = v1(range4);
    v2_1 = v2(range4);
    p = w1ps(:,range6)*wt;
    ic2(range4) = ic2(range4)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    p = w2ps(:,range6)*wt;
    ic3(range4) = ic3(range4)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    range5 = lds.phases;         
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        sjac  = cjac(lds.func,lds.Jacobian,xt,pt,lds.ActiveParams);
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        sysjac(range5,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',sjac);

        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shv2(:,d1) = (wtk.*hess(:,kr2,d1))*v2_1;
        end
        fxhessv1v1(range2)=wtk.*shv1(:,kr2)*v1_1;
        fxhessv1v2(range2)=wtk.*shv1(:,kr2)*v2_1;
        fxhessv2v2(range2)=wtk.*shv2(:,kr2)*v2_1;
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    fxjac(range3,range4)=sysjac;
    range1 = range1+lds.ncol;
    range3 = range3 + lds.ncol_coord;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end

%computation of a20
a20 = 1/2*phi'*(fxhessv1v1'+2*fxjac*v1);

%computation of alpha20
%alpha20 = 1/2*w1'*(fxhessv1v1'+2*fxjac*v1)+1;
alpha20 = 0;

%computation of h20
jace(lds.ncoords+1,lds.coords) = ic2;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v1'-2*T*a20*v1C-2*T*alpha20*f1+2*T*fxjac*v1+2*T*f1; rl(lds.ncoords+1) = 0;
h20 = jace\rl';
h20 = h20(lds.coords);

%computation of b11
b11 = w2'*(fxhessv1v2'+fxjac*v2);

%computation of h11
Jjac(lds.ncoords+1,lds.coords) = ic3;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v2'-T*b11*v2C+T*fxjac*v2; rl(lds.ncoords+1) = 0;
h11 = Jjac\rl';
h11 = h11(lds.coords);

%computation of a021
a021 = 1/2*phi'*fxhessv2v2';
a02 = a021/T;

%computation of alpha021
%alpha021 = 1/2*w1'*fxhessv2v2';
alpha021 = 0;

rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv2v2'-2*a021*T*v1C-2*alpha021*T*f1; rl(lds.ncoords+1) = 0;
h02 = jace\rl';
h02 = h02(lds.coords);

% function
range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
h20ps = reshape(h20,lds.nphase,lds.tps);
h02ps = reshape(h02,lds.nphase,lds.tps);
h11ps = reshape(h11,lds.nphase,lds.tps);
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    h20Cps(:,range6) = h20ps(:,range1)*lds.wt;
    h02Cps(:,range6) = h02ps(:,range1)*lds.wt;
    h11Cps(:,range6) = h11ps(:,range1)*lds.wt;
    v1_1 = v1(range4);
    v11_1 = conj(v1_1);
    v2_1 = v2(range4);
    h02_1 = h02(range4);
    h11_1 = h11(range4);
    h20_1 = h20(range4);
    range5 = lds.phases;         
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        stenv1v1 = zeros(lds.nphase,lds.nphase);
        stenv1v2 = zeros(lds.nphase,lds.nphase);
        stenv2v2 = zeros(lds.nphase,lds.nphase);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        for d1=lds.phases
            shh02(:,d1) = (wtk.*hess(:,kr2,d1))*h02_1;
            shh20(:,d1) = (wtk.*hess(:,kr2,d1))*h20_1;
            shh11(:,d1) = (wtk.*hess(:,kr2,d1))*h11_1;
            for d2=lds.phases
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
                stensv2(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v2_1;
            end    
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
            stenv1v2(:,d1)=stenv1v2(:,d1)+wtk.*stensv1(:,kr2,d1)*v2_1;
            stenv2v2(:,d1)=stenv2v2(:,d1)+wtk.*stensv2(:,kr2,d1)*v2_1;
        end
        fxhessh02v1(range2)=wtk.*shh02(:,kr2)*v1_1;
        fxhessh02v2(range2)=wtk.*shh02(:,kr2)*v2_1;
        fxhessh20v1(range2)=wtk.*shh20(:,kr2)*v1_1;
        fxhessh20v2(range2)=wtk.*shh20(:,kr2)*v2_1;        
        fxhessh11v11(range2)=wtk.*shh11(:,kr2)*v11_1;
        fxhessh11v2(range2)=wtk.*shh11(:,kr2)*v2_1;
        fxtensv1v1v1(range2)=wtk.*stenv1v1(:,kr2)*v1_1;
        fxtensv1v1v2(range2)=wtk.*stenv1v1(:,kr2)*v2_1;
        fxtensv1v2v2(range2)=wtk.*stenv1v2(:,kr2)*v2_1;
        fxtensv2v2v2(range2)=wtk.*stenv2v2(:,kr2)*v2_1;
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end
h20C = reshape(h20Cps,lds.ncoords-lds.nphase,1);
h02C = reshape(h02Cps,lds.ncoords-lds.nphase,1);
h11C = reshape(h11Cps,lds.ncoords-lds.nphase,1);

%computation of a30
a30 = 1/6*phi'*(fxtensv1v1v1'+3*fxhessh20v1'-6*a20*h20C+3*(fxjac*h20+fxhessv1v1')+6*(1-alpha20)*fxjac*v1)-a20;

%computation of b21
b21 = 1/2*w2'*(fxtensv1v1v2'+fxhessh20v2'+2*fxhessh11v11'-2*a20*h11C-2*b11*h11C+2*(fxjac*h11+fxhessv1v2')+2*(1-alpha20)*fxjac*v2)-b11;

%computation of a12
a12 = 1/(2*T)*phi'*(fxtensv1v2v2'+fxhessh02v1'+2*fxhessh11v2'-2*b11*h02C-2*a021*h20C+fxjac*h02+fxhessv2v2'-2*alpha021*fxjac*v1)-a021/T;

%computation of b03
b03 = 1/(6*T)*w2'*(fxtensv2v2v2'+3*fxhessh02v2'-6*a021*h11C-6*alpha021*fxjac*v2);

a1 = -a20/b11;
b1 = -a021/(T*b11);
%lns = -2*a20^2*b03-3*a021/T*a30*b11 + a20*(a12*b11 + 2*b03*b11 + 2*a021/T*b21);
cns = -2*a20*b21*a02+6*b03*a20^2+(-a02*b21-6*a20*a02+2*a20*b03-3*a02*a30-a12*a20)*b11+b11^2*(a12-a02);
result = [b11 a20 a021/T cns];

% ---------------------------------------------------------------
function [x,p,T] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global lds

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+1:lds.PeriodIdx+2);
x = x0(lds.coords);
T = x0(lds.PeriodIdx);
