function result=nf_PDNS(x)
%
% calculates PDNS normal form coefficient
%
global lds

higherorderderivatives = 0;

lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);

[xt,p,T] = rearr(x);
pt = num2cell(p);
ups = reshape(xt,lds.nphase,lds.tps);

d = lds.multipliers;
smallest_sum = Inf;
    for jk=1:lds.nphase-1
        if abs(imag(d(jk))) > 1e-3        
            [val,idx] = min(abs(d(jk+1:lds.nphase)*d(jk)-1));
            if val < smallest_sum
                idx2 = jk+idx;
                smallest_sum = val;
            end
        end
    end
theta = abs(angle(d(idx2)));

%computation main part matrices
[Jjac,Jjac2] = BVP_jac('bordBVP_R2_f','BVP_LC1_jac_bc',xt,p,T,theta,1,1);
Jjac3 = Jjac2;
Jjac3(lds.ncoords-lds.nphase+1:lds.ncoords,[lds.phases lds.ncoords-lds.nphase+1:lds.ncoords]) = bordBVP_R2_bc1;

%computation of v1
b = []; b(lds.ncoords+1)=1; b=b';
jace=[Jjac3,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

Jjac3=[Jjac3 p;q' 0];
[LJ,UJ] = lu(Jjac3);
vext = UJ \(LJ\b);

%rescale vext
v1 = vext(1:lds.ncoords);
vps = reshape(v1,lds.nphase,lds.tps);
ficd = dot(vps,vps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v1= v1/sqrt(ic);

%computation of v2
jace=[Jjac,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

Jjac=[Jjac p;q' 0];
[LJ,UJ] = lu(Jjac);
vext = UJ \(LJ\b);

%rescale vext
v2 = vext(1:lds.ncoords);
vps = reshape(v2,lds.nphase,lds.tps);
ficd = dot(vps,vps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v2= v2/sqrt(ic);

jace=[Jjac2,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

Jjac2=[Jjac2 p;q' 0];

%computation of phi*
phi = Jjac2'\b;
phi = phi(1:lds.ncoords-lds.nphase)';

%computation of v1C,v2C,f1
v1ps = reshape(v1,lds.nphase,lds.tps);
v2ps = reshape(v2,lds.nphase,lds.tps);
range1 = lds.cols_p1;
range2 = 1 : lds.ncol;
range3 = lds.phases;
for j=lds.tsts
    v1Cps(:,range2) = v1ps(:,range1)*lds.wt;
    v2Cps(:,range2) = v2ps(:,range1)*lds.wt;
    xpp = ups(:,range1)*lds.wt;
    for c=lds.cols                    
        xt = xpp(:,c);
        f1(range3) = feval(lds.func, 0,xt, pt{:});
        range3 = range3 + lds.nphase;
    end   
    range1 = range1 + lds.ncol;
    range2 = range2 + lds.ncol;
end
v1C = reshape(v1Cps,lds.ncoords-lds.nphase,1);
v2C = reshape(v2Cps,lds.ncoords-lds.nphase,1);
f1 = f1';

%rescale phi*
phi = phi/(phi*f1);

%computation of v1*
w1 = Jjac3'\b;
w1 = w1(1:lds.ncoords-lds.nphase)';
w1 = w1/(w1*v1C);

%computation of v2*
w2 = Jjac'\b;
w2 = w2(1:lds.ncoords-lds.nphase)';
w2 = w2/(w2*v2C);

% function
range1 = lds.cols_p1;
range2 = lds.phases;
range3 = lds.col_coords;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
phips = reshape(phi,lds.nphase,lds.tps-1);
w1ps = reshape(w1,lds.nphase,lds.tps-1);
w2ps = reshape(w2,lds.nphase,lds.tps-1);
w22ps = conj(w2ps);
icphi = zeros(1,lds.ncoords);
icw1 = zeros(1,lds.ncoords);
icw22 = zeros(1,lds.ncoords);
wt = lds.wt';
for j=lds.tsts
    % value of polynomial on each collocation point
    xp(:,range6) = ups(:,range1)*lds.wt;
    xpp = xp(:,range6);
    v1_1 = v1(range4);
    v2_1 = v2(range4);
    v22_1 = conj(v2_1);
    p = phips(:,range6)*wt;
    icphi(range4) = icphi(range4)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    p = w1ps(:,range6)*wt;
    icw1(range4) = icw1(range4)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    p = w22ps(:,range6)*wt;
    icw22(range4) = icw22(range4)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    range5 = lds.phases;         
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xpp(:,c); 
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
        fxhessv2v22(range2)=wtk.*shv2(:,kr2)*v22_1;
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    fxjac(range3,range4)=sysjac;
    range1 = range1+lds.ncol;
    range3 = range3 + lds.ncol_coord;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;    
end

%computation of alpha2001
alpha2001 = 1/2*phi*fxhessv1v1';
alpha200 = alpha2001/T;

%computation of h200
Jjac2(lds.ncoords+1,lds.coords) = icphi;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v1'-2*T*alpha2001*f1; rl(lds.ncoords+1) = 0;
diffh200 = rl(1:lds.ncoords-lds.nphase);
h200 = Jjac2\rl.';
h200 = h200(lds.coords);

%computation of h020
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
range3 = lds.cols;
twitheta = 2*1i*theta*lds.pwwt;
J2thetamin = Jjac;
for jk  = lds.tsts
  % evaluate part of Jacobian matrix
  J2thetamin(range1,range2) = bordBVP_R2_f(lds.func,xp(:,range3),pt,T,jk)+twitheta;  
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
  range3 = range3 + lds.ncol;
end
% remove borders
J2thetamin(:,end)=[];J2thetamin(end,:)=[];
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv2v2.'; rl(lds.ncoords) = 0;
diffh020 = rl(1:lds.ncoords-lds.nphase);
h020 = J2thetamin\rl.';

%computation of h110
Jthetaplus = Jjac;
Jthetaplus(lds.ncoords-lds.nphase+1:lds.ncoords,[lds.phases lds.ncoords-lds.nphase+1:lds.ncoords]) = bordBVP_R2_bc1;
Jthetaplus(:,end) = [];Jthetaplus(end,:) = [];
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v2.'; rl(lds.ncoords) = 0;
diffh110 = rl(1:lds.ncoords-lds.nphase);
h110 = Jthetaplus\rl.';

%computation of alpha0111
alpha0111 = phi*real(fxhessv2v22)';
alpha011 = alpha0111/T;

%computation of h011
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxhessv2v22)'-T*alpha0111*f1; rl(lds.ncoords+1) = 0;
diffh011 = rl(1:lds.ncoords-lds.nphase);
h011 = Jjac2\rl.';
h011 = h011(lds.coords);

% function
range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    xpp(:,range6) = xp;
    v1_1 = v1(range4);
    v2_1 = v2(range4);
    v22_1 = conj(v2_1);
    h110_1 = h110(range4);
    h101_1 = conj(h110_1);
    h200_1 = h200(range4);
    h011_1 = h011(range4);
    h020_1 = h020(range4);
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
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shv2(:,d1) = (wtk.*hess(:,kr2,d1))*v2_1;
            shv22(:,d1) = (wtk.*hess(:,kr2,d1))*v22_1;
            for d2=lds.phases
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
                stensv2(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v2_1;
            end    
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
            stenv1v2(:,d1)=stenv1v2(:,d1)+wtk.*stensv1(:,kr2,d1)*v2_1;
            stenv2v2(:,d1)=stenv2v2(:,d1)+wtk.*stensv2(:,kr2,d1)*v2_1;
        end
        fxhessv1h200(range2)=wtk.*shv1(:,kr2)*h200_1;
        fxhessv1h110(range2)=wtk.*shv1(:,kr2)*h110_1;
        fxhessv2h200(range2)=wtk.*shv2(:,kr2)*h200_1;
        fxhessv1h011(range2)=wtk.*shv1(:,kr2)*h011_1;
        fxhessv2h011(range2)=wtk.*shv2(:,kr2)*h011_1;
        fxhessv1h020(range2)=wtk.*shv1(:,kr2)*h020_1;
        fxhessv22h020(range2)=wtk.*shv22(:,kr2)*h020_1;
        fxhessv2h110(range2)=wtk.*shv2(:,kr2)*h110_1;
        fxhessv2h101(range2)=wtk.*shv2(:,kr2)*h101_1;
        fxhessv2h020(range2)=wtk.*shv2(:,kr2)*h020_1;        
        fxtensv1v1v1(range2)=wtk.*stenv1v1(:,kr2)*v1_1;
        fxtensv1v1v2(range2)=wtk.*stenv1v1(:,kr2)*v2_1;
        fxtensv1v2v2(range2)=wtk.*stenv1v2(:,kr2)*v2_1;
        fxtensv1v2v22(range2)=wtk.*stenv1v2(:,kr2)*v22_1;
        fxtensv2v2v2(range2)=wtk.*stenv2v2(:,kr2)*v2_1;
        fxtensv2v2v22(range2)=wtk.*stenv2v2(:,kr2)*v22_1;
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
end

%computation of a3001
a3001 = 1/6*w1*(fxtensv1v1v1'+3*fxhessv1h200'-6*alpha2001*fxjac*v1);
a300 = a3001/T;

%computation of h300
Jjac3(lds.ncoords+1,lds.coords) = icw1;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v1v1'+3*T*fxhessv1h200'-6*T*alpha2001*fxjac*v1-6*T*a3001*v1C; rl(lds.ncoords+1) = 0;
diffh300 = rl(1:lds.ncoords-lds.nphase);
h300 = Jjac3\rl.';
h300 = h300(lds.coords);

%computation of h030
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
range3 = lds.cols;
threeitheta = 3*1i*theta*lds.pwwt;
J3thetamin = Jjac;
for jk  = lds.tsts
  % evaluate part of Jacobian matrix 
  J3thetamin(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+threeitheta;  
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
  range3 = range3 + lds.ncol;
end

% remove borders
J3thetamin(:,end)=[];J3thetamin(end,:)=[];
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv2v2v2.'+3*T*fxhessv2h020.'; rl(lds.ncoords) = 0;
h030 = J3thetamin\rl.';

%computation of b2101
b2101 = 1/2*w2*(fxtensv1v1v2.'+fxhessv2h200.'+2*fxhessv1h110.'-2*alpha2001*fxjac*v2)+i*theta/T*alpha2001;
b210 = b2101/T;

%computation of h210
Jjac(lds.ncoords+1,lds.coords) = conj(icw22);
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v1v2.'+T*fxhessv2h200.'+2*T*fxhessv1h110.'-2*T*alpha2001*fxjac*v2...
    +2*alpha2001*i*theta*v2C-2*T*b2101*v2C; rl(lds.ncoords+1) = 0;
diffh210 = rl(1:lds.ncoords-lds.nphase);
h210 = Jjac\rl.';
h210 = h210(lds.coords);

%computation of h120
J2thetaplus = J2thetamin;
J2thetaplus(lds.ncoords-lds.nphase+1:lds.ncoords,[lds.phases lds.ncoords-lds.nphase+1:lds.ncoords]) = bordBVP_R2_bc1;
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v2v2.'+T*fxhessv1h020.'+2*T*fxhessv2h110.'; rl(lds.ncoords) = 0;
h120 = J2thetaplus\rl.';

%computation of b0211
b0211 = 1/2*w2*(fxtensv2v2v22.'+fxhessv22h020.'+2*fxhessv2h011.'-2*alpha0111*fxjac*v2)+alpha0111*i*theta/T;
b021 = b0211/T;

%computation of a1111
a1111 = w1*(real(fxtensv1v2v22)'+fxhessv1h011'+2*real(fxhessv2h101)'-alpha0111*fxjac*v1);
a111 = a1111/T;

%computation of h021
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv2v2v22.'+T*fxhessv22h020.'+2*T*fxhessv2h011.'-2*T*alpha0111*fxjac*v2...
    +2*alpha0111*i*theta*v2C-2*T*b0211*v2C; rl(lds.ncoords+1) = 0;
diffh021 = rl(1:lds.ncoords-lds.nphase);
h021 = Jjac\rl.';
h021 = h021(lds.coords);

%computation of h111
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxtensv1v2v22)'+T*fxhessv1h011'+2*T*real(fxhessv2h101')-T*alpha0111*fxjac*v1...
    -a1111*T*v1C; rl(lds.ncoords+1) = 0;
diffh111 = rl(1:lds.ncoords-lds.nphase);
h111 = Jjac3\rl.';
h111 = h111(lds.coords);

if higherorderderivatives == 1
range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
h200ps = reshape(h200,lds.nphase,lds.tps);
h020ps = reshape(h020,lds.nphase,lds.tps);
h210ps = reshape(h210,lds.nphase,lds.tps);
h021ps = reshape(h021,lds.nphase,lds.tps);
h110ps = reshape(h110,lds.nphase,lds.tps);
h011ps = reshape(h011,lds.nphase,lds.tps);
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    h200Cps(:,range6) = h200ps(:,range1)*lds.wt;
    h020Cps(:,range6) = h020ps(:,range1)*lds.wt;
    h210Cps(:,range6) = h210ps(:,range1)*lds.wt;
    h021Cps(:,range6) = h021ps(:,range1)*lds.wt;
    h110Cps(:,range6) = h110ps(:,range1)*lds.wt;
    h011Cps(:,range6) = h011ps(:,range1)*lds.wt;
    v1_1 = v1(range4);
    v2_1 = v2(range4);
    v22_1 = conj(v2_1);    
    h011_1 = h011(range4);
    h021_1 = h021(range4);
    h012_1 = conj(h021_1);
    h020_1 = h020(range4);
    h002_1 = conj(h020_1);
    h030_1 = h030(range4);   
    h110_1 = h110(range4);
    h101_1 = conj(h110_1);
    h111_1 = h111(range4);
    h120_1 = h120(range4);
    h200_1 = h200(range4);     
    h210_1 = h210(range4);
    h201_1 = conj(h210_1);
    h300_1 = h300(range4);
    
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        stenv1v1 = zeros(lds.nphase,lds.nphase);
        stenv1v2 = zeros(lds.nphase,lds.nphase);
        stenv1v22 = zeros(lds.nphase,lds.nphase);
        stenv2v2 = zeros(lds.nphase,lds.nphase);
        stenv2v22 = zeros(lds.nphase,lds.nphase);        
        stens4v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v1v2 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v2v2 = zeros(lds.nphase,lds.nphase,lds.nphase);        
        stens4v1v1v1 = zeros(lds.nphase,lds.nphase);
        stens4v1v1v2 = zeros(lds.nphase,lds.nphase);
        stens4v1v2v2 = zeros(lds.nphase,lds.nphase);
        stens4v2v2v2 = zeros(lds.nphase,lds.nphase);
        stens4v2v2v22 = zeros(lds.nphase,lds.nphase);        
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        tens4 = ctens4(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,xt,pt,lds.ActiveParams);

        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shv2(:,d1) = (wtk.*hess(:,kr2,d1))*v2_1;
            shv22(:,d1) = (wtk.*hess(:,kr2,d1))*v22_1;
            shh020(:,d1) = (wtk.*hess(:,kr2,d1))*h020_1;
            shh200(:,d1) = (wtk.*hess(:,kr2,d1))*h200_1;
            shh101(:,d1) = (wtk.*hess(:,kr2,d1))*h101_1;
            shh110(:,d1) = (wtk.*hess(:,kr2,d1))*h110_1;
            shh011(:,d1) = (wtk.*hess(:,kr2,d1))*h011_1;
            for d2=lds.phases
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
                stensv2(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v2_1;
                for d3=lds.phases
                    stens4v1(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v1_1;
                    stens4v2(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v2_1;
                end
                stens4v1v1(:,d2,d1) = stens4v1v1(:,d2,d1)+(wtk.*stens4v1(:,kr2,d2,d1))*v1_1;
                stens4v1v2(:,d2,d1) = stens4v1v2(:,d2,d1)+(wtk.*stens4v1(:,kr2,d2,d1))*v2_1;
                stens4v2v2(:,d2,d1) = stens4v2v2(:,d2,d1)+(wtk.*stens4v2(:,kr2,d2,d1))*v2_1;
            end  
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
            stenv1v2(:,d1)=stenv1v2(:,d1)+wtk.*stensv1(:,kr2,d1)*v2_1;
            stenv1v22(:,d1)=stenv1v22(:,d1)+wtk.*stensv1(:,kr2,d1)*v22_1;
            stenv2v2(:,d1)=stenv2v2(:,d1)+wtk.*stensv2(:,kr2,d1)*v2_1;
            stenv2v22(:,d1)=stenv2v22(:,d1)+wtk.*stensv2(:,kr2,d1)*v22_1;
            stens4v1v1v1(:,d1) = stens4v1v1v1(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v1_1;
            stens4v1v1v2(:,d1) = stens4v1v1v2(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v2_1;
            stens4v1v2v2(:,d1) = stens4v1v2v2(:,d1)+(wtk.*stens4v1v2(:,kr2,d1))*v2_1;
            stens4v2v2v2(:,d1) = stens4v2v2v2(:,d1)+(wtk.*stens4v2v2(:,kr2,d1))*v2_1;
            stens4v2v2v22(:,d1) = stens4v2v2v22(:,d1)+(wtk.*stens4v2v2(:,kr2,d1))*v22_1;
        end
        fxhessv1h300(range2)=wtk.*shv1(:,kr2)*h300_1;
        fxhessv1h030(range2)=wtk.*shv1(:,kr2)*h030_1;
        fxhessv1h210(range2)=wtk.*shv1(:,kr2)*h210_1;
        fxhessv1h120(range2)=wtk.*shv1(:,kr2)*h120_1;
        fxhessv1h021(range2)=wtk.*shv1(:,kr2)*h021_1;
        fxhessv1h111(range2)=wtk.*shv1(:,kr2)*h111_1;
        fxhessv2h030(range2)=wtk.*shv2(:,kr2)*h030_1;
        fxhessv2h300(range2)=wtk.*shv2(:,kr2)*h300_1;
        fxhessv2h120(range2)=wtk.*shv2(:,kr2)*h120_1;
        fxhessv2h021(range2)=wtk.*shv2(:,kr2)*h021_1;
        fxhessv2h012(range2)=wtk.*shv2(:,kr2)*h012_1;
        fxhessv2h201(range2)=wtk.*shv2(:,kr2)*h201_1;
        fxhessv2h210(range2)=wtk.*shv2(:,kr2)*h210_1;
        fxhessv2h111(range2)=wtk.*shv2(:,kr2)*h111_1;
        fxhessv22h030(range2)=wtk.*shv22(:,kr2)*h030_1;
        fxhessv22h120(range2)=wtk.*shv22(:,kr2)*h120_1;
        fxhessh200h200(range2)=wtk.*shh200(:,kr2)*h200_1;
        fxhessh200h110(range2)=wtk.*shh200(:,kr2)*h110_1;
        fxhessh200h011(range2)=wtk.*shh200(:,kr2)*h011_1;
        fxhessh200h020(range2)=wtk.*shh200(:,kr2)*h020_1;
        fxhessh020h020(range2)=wtk.*shh020(:,kr2)*h020_1;
        fxhessh020h002(range2)=wtk.*shh020(:,kr2)*h002_1;
        fxhessh020h110(range2)=wtk.*shh020(:,kr2)*h110_1;
        fxhessh020h101(range2)=wtk.*shh020(:,kr2)*h101_1;
        fxhessh020h011(range2)=wtk.*shh020(:,kr2)*h011_1;
        fxhessh101h110(range2)=wtk.*shh101(:,kr2)*h110_1;
        fxhessh110h110(range2)=wtk.*shh110(:,kr2)*h110_1;
        fxhessh011h110(range2)=wtk.*shh011(:,kr2)*h110_1;
        fxhessh011h011(range2)=wtk.*shh011(:,kr2)*h011_1;
        fxtensv1v1h200(range2)=wtk.*stenv1v1(:,kr2)*h200_1;
        fxtensv1v1h020(range2)=wtk.*stenv1v1(:,kr2)*h020_1;
        fxtensv1v1h110(range2)=wtk.*stenv1v1(:,kr2)*h110_1;
        fxtensv1v1h011(range2)=wtk.*stenv1v1(:,kr2)*h011_1;
        fxtensv1v2h200(range2)=wtk.*stenv1v2(:,kr2)*h200_1;
        fxtensv1v2h020(range2)=wtk.*stenv1v2(:,kr2)*h020_1;
        fxtensv1v2h101(range2)=wtk.*stenv1v2(:,kr2)*h101_1;
        fxtensv1v2h110(range2)=wtk.*stenv1v2(:,kr2)*h110_1;
        fxtensv1v2h011(range2)=wtk.*stenv1v2(:,kr2)*h011_1;
        fxtensv1v22h020(range2)=wtk.*stenv1v22(:,kr2)*h020_1;
        fxtensv2v2h020(range2)=wtk.*stenv2v2(:,kr2)*h020_1;
        fxtensv2v2h002(range2)=wtk.*stenv2v2(:,kr2)*h002_1;
        fxtensv2v2h110(range2)=wtk.*stenv2v2(:,kr2)*h110_1;
        fxtensv2v2h200(range2)=wtk.*stenv2v2(:,kr2)*h200_1;
        fxtensv2v2h101(range2)=wtk.*stenv2v2(:,kr2)*h101_1;
        fxtensv2v2h011(range2)=wtk.*stenv2v2(:,kr2)*h011_1;
        fxtensv2v22h020(range2)=wtk.*stenv2v22(:,kr2)*h020_1;
        fxtensv2v22h200(range2)=wtk.*stenv2v22(:,kr2)*h200_1;
        fxtensv2v22h110(range2)=wtk.*stenv2v22(:,kr2)*h110_1;
        fxtensv2v22h011(range2)=wtk.*stenv2v22(:,kr2)*h011_1;
        fxtens4v1v1v1v1(range2)=wtk.*stens4v1v1v1(:,kr2)*v1_1;
        fxtens4v1v1v1v2(range2)=wtk.*stens4v1v1v1(:,kr2)*v2_1;
        fxtens4v1v1v2v2(range2)=wtk.*stens4v1v1v2(:,kr2)*v2_1;
        fxtens4v1v1v2v22(range2)=wtk.*stens4v1v1v2(:,kr2)*v22_1;
        fxtens4v1v2v2v2(range2)=wtk.*stens4v1v2v2(:,kr2)*v2_1;        
        fxtens4v1v2v2v22(range2)=wtk.*stens4v1v2v2(:,kr2)*v22_1;        
        fxtens4v2v2v2v2(range2)=wtk.*stens4v2v2v2(:,kr2)*v2_1;        
        fxtens4v2v2v22v22(range2)=wtk.*stens4v2v2v22(:,kr2)*v22_1;        
        fxtens4v2v2v2v22(range2)=wtk.*stens4v2v2v2(:,kr2)*v22_1;        
        range2 = range2+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end
h200C = reshape(h200Cps,lds.ncoords-lds.nphase,1);
h020C = reshape(h020Cps,lds.ncoords-lds.nphase,1);
h210C = reshape(h210Cps,lds.ncoords-lds.nphase,1);
h021C = reshape(h021Cps,lds.ncoords-lds.nphase,1);
h110C = reshape(h110Cps,lds.ncoords-lds.nphase,1);
h011C = reshape(h011Cps,lds.ncoords-lds.nphase,1);

%computation of alpha4001
diffh200 = T*fxjac*h200+diffh200';
alpha4001 = 1/24*phi*(fxtens4v1v1v1v1'+6*fxtensv1v1h200'+3*fxhessh200h200'+4*fxhessv1h300'-12*alpha2001*diffh200);
alpha400 = alpha4001/T^2;

%computation of h400
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v1v1'+6*T*fxtensv1v1h200'+3*T*fxhessh200h200'+4*T*fxhessv1h300'...
    -12*T*alpha2001*diffh200-24*T*alpha4001*f1-24*T*a3001*h200C; rl(lds.ncoords+1) = 0;
h400 = Jjac2\rl.';
h400 = h400(lds.coords);

%computation of h040
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
range3 = lds.cols;
fouritheta = 4*1i*theta*lds.pwwt;
J4thetamin = Jjac;
for jk  = lds.tsts
  % evaluate part of Jacobian matrix
  J4thetamin(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+fouritheta;  
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
  range3 = range3 + lds.ncol;
end
% remove borders
J4thetamin(:,end)=[];J4thetamin(end,:)=[];
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v2v2v2v2.'+6*T*fxtensv2v2h020.'+4*T*fxhessv2h030.'+3*T*fxhessh020h020.'; rl(lds.ncoords) = 0;
h040 = J4thetamin\rl.';

%computation of h310
diffh110 = T*fxjac*h110-i*theta*h110C+diffh110.';
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v1v2.'+3*T*fxtensv1v1h110.'+3*T*fxtensv1v2h200.'+T*fxhessv2h300.'...
    +3*T*fxhessv1h210.'+3*T*fxhessh200h110.'-6*T*alpha2001*diffh110-6*T*a3001*h110C-6*T*b2101*h110C; rl(lds.ncoords) = 0;
h310 = Jthetaplus\rl.';

%computation of h130
J3thetaplus = J3thetamin;
J3thetaplus(lds.ncoords-lds.nphase+1:lds.ncoords,[lds.phases lds.ncoords-lds.nphase+1:lds.ncoords]) = bordBVP_R2_bc1;
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v2v2v2.'+3*T*fxtensv2v2h110.'+3*T*fxtensv1v2h020.'+T*fxhessv1h030.'...
    +3*T*fxhessh020h110.'+3*T*fxhessv2h120.'; rl(lds.ncoords) = 0;
h130 = J3thetaplus\rl.';

%computation of h031
diffh020 = T*fxjac*h020-2*i*theta*h020C+diffh020.';
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v2v2v2v22.'+3*T*fxtensv2v2h011.'+3*T*fxtensv2v22h020.'+T*fxhessv22h030.'...
    +3*T*fxhessh020h011.'+3*T*fxhessv2h021.'-3*T*alpha0111*diffh020-6*T*b0211*h020C; rl(lds.ncoords) = 0;
h031 = J2thetamin\rl.';

%computation of alpha2111
diffh011 = T*fxjac*h011+diffh011.';
alpha2111 = 1/2*phi*(real(fxtens4v1v1v2v22)'+fxtensv1v1h011'+real(fxtensv2v22h200)'+4*real(fxtensv1v2h101')...
    +2*real(fxhessv2h201')+fxhessh200h011'+2*real(fxhessh101h110)'+2*fxhessv1h111'-alpha0111*diffh200-2*alpha2001*diffh011);
alpha211 = alpha2111/T^2;

%computation of h211
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxtens4v1v1v2v22)'+T*fxtensv1v1h011'+T*real(fxtensv2v22h200)'...
    +4*T*real(fxtensv1v2h101)'+2*T*real(fxhessv2h201)'+T*fxhessh200h011'+2*T*real(fxhessh101h110)'...
    +2*T*fxhessv1h111'-T*alpha0111*diffh200-2*T*alpha2001*diffh011-2*T*alpha2111*f1-2*T*a1111*h200C-4*T*real(b2101)*h011C; rl(lds.ncoords+1) = 0;
h211 = Jjac2\rl.';
h211 = h211(lds.coords);

%computation of h121
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v2v2v22.'+T*fxtensv1v22h020.'+2*T*fxtensv1v2h011.'+T*fxtensv2v2h101.'...
    +2*T*fxtensv2v22h110.'+T*fxhessv1h021.'+T*fxhessh020h101.'+2*T*fxhessh011h110.'+2*T*fxhessv2h111.'...
    +T*fxhessv22h120.'-2*T*alpha0111*diffh110-2*T*b0211*h110C-2*T*a1111*h110C; rl(lds.ncoords) = 0;
h121 = Jthetaplus\rl.';

%computation of h220
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v2v2.'+T*fxtensv2v2h200.'+4*T*fxtensv1v2h110.'+T*fxtensv1v1h020.'...
    +T*fxhessh200h020.'+2*T*fxhessv2h210.'+2*T*fxhessh110h110.'+2*T*fxhessv1h120.'-2*T*alpha2001*diffh020-4*T*b2101*h020C; rl(lds.ncoords) = 0;
h220 = J2thetamin\rl.';

%computation of alpha0221
alpha0221 = 1/4*phi*(real(fxtens4v2v2v22v22)'+4*real(fxtensv2v22h011)'+2*real(fxtensv2v2h002)'...
    +real(fxhessh020h002)'+2*fxhessh011h011'+4*real(fxhessv2h012)'-4*alpha0111*diffh011);
alpha022 = alpha0221/T^2;

%computation of h022
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxtens4v2v2v22v22)'+T*4*real(fxtensv2v22h011)'+T*2*real(fxtensv2v2h002)'...
    +T*real(fxhessh020h002)'+T*2*fxhessh011h011'+4*T*real(fxhessv2h012)'-4*T*alpha0111*diffh011-4*T*alpha0221*f1...
    -8*real(b0211)*h011C*T; rl(lds.ncoords+1) = 0;
h022 = Jjac2\rl.';
h022 = h022(lds.coords);

range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    v1_1 = v1(range4);
    v2_1 = v2(range4);
    v22_1 = conj(v2_1);        
    h011_1 = h011(range4);    
    h020_1 = h020(range4);
    h002_1 = conj(h020_1);
    h021_1 = h021(range4);
    h012_1 = conj(h021_1);
    h022_1 = h022(range4);
    h030_1 = h030(range4);
    h031_1 = h031(range4);        
    h110_1 = h110(range4);
    h101_1 = conj(h110_1);
    h111_1 = h111(range4);    
    h120_1 = h120(range4);    
    h102_1 = conj(h120_1);
    h121_1 = h121(range4);
    h112_1 = conj(h121_1);
    h200_1 = h200(range4);    
    h210_1 = h210(range4);
    h201_1 = conj(h210_1);
    h211_1 = h211(range4);
    h220_1 = h220(range4);
    h300_1 = h300(range4);    
    h310_1 = h310(range4);
    h301_1 = conj(h310_1);
    h400_1 = h400(range4);
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        stens5v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase,lds.nphase);
        stens5v1v2 = zeros(lds.nphase,lds.nphase,lds.nphase,lds.nphase);
        stens5v2v2 = zeros(lds.nphase,lds.nphase,lds.nphase,lds.nphase);        
        sten4v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase);
        sten4v1v2 = zeros(lds.nphase,lds.nphase,lds.nphase);
        sten4v2v2 = zeros(lds.nphase,lds.nphase,lds.nphase);
        sten4v2v22 = zeros(lds.nphase,lds.nphase,lds.nphase);        
        stens5v1v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase); 
        stens5v1v1v2 = zeros(lds.nphase,lds.nphase,lds.nphase); 
        stens5v1v2v2 = zeros(lds.nphase,lds.nphase,lds.nphase); 
        stens5v2v2v2 = zeros(lds.nphase,lds.nphase,lds.nphase); 
        
        stenv1v1 = zeros(lds.nphase,lds.nphase);
        stenv1v2 = zeros(lds.nphase,lds.nphase);
        stenv1v22 = zeros(lds.nphase,lds.nphase);
        stenv1h200 = zeros(lds.nphase,lds.nphase); 
        stenv1h020 = zeros(lds.nphase,lds.nphase); 
        stenv1h101 = zeros(lds.nphase,lds.nphase); 
        stenv1h011 = zeros(lds.nphase,lds.nphase); 
        stenv1h110 = zeros(lds.nphase,lds.nphase); 
        stenv2v2 = zeros(lds.nphase,lds.nphase); 
        stenv22v22 = zeros(lds.nphase,lds.nphase); 
        stenv2h200 = zeros(lds.nphase,lds.nphase); 
        stenv2h101 = zeros(lds.nphase,lds.nphase); 
        stenv2h002 = zeros(lds.nphase,lds.nphase); 
        stenv2h011 = zeros(lds.nphase,lds.nphase); 
        stenv22h110 = zeros(lds.nphase,lds.nphase); 
        stenv22h101 = zeros(lds.nphase,lds.nphase); 
        stenv22h200 = zeros(lds.nphase,lds.nphase); 
        stenv22h020 = zeros(lds.nphase,lds.nphase); 
            
        sten4v1v1v1 = zeros(lds.nphase,lds.nphase); 
        sten4v1v1v2 = zeros(lds.nphase,lds.nphase); 
        sten4v1v1v22 = zeros(lds.nphase,lds.nphase); 
        sten4v1v2v2 = zeros(lds.nphase,lds.nphase); 
        sten4v1v2v22 = zeros(lds.nphase,lds.nphase); 
        sten4v2v2v2 = zeros(lds.nphase,lds.nphase); 
        sten4v2v2v22 = zeros(lds.nphase,lds.nphase); 
        sten4v2v22v22 = zeros(lds.nphase,lds.nphase); 
        stens5v1v1v1v1 = zeros(lds.nphase,lds.nphase); 
        stens5v1v1v1v2 = zeros(lds.nphase,lds.nphase); 
        stens5v1v1v2v2 = zeros(lds.nphase,lds.nphase); 
        stens5v1v2v2v22 = zeros(lds.nphase,lds.nphase); 
        stens5v2v2v2v22 = zeros(lds.nphase,lds.nphase);              
        
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        tens4 = ctens4(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,xt,pt,lds.ActiveParams);
        tens5 = ctens5(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,lds.Der5,xt,pt,lds.ActiveParams);
        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shv2(:,d1) = (wtk.*hess(:,kr2,d1))*v2_1;
            shv22(:,d1) = (wtk.*hess(:,kr2,d1))*v22_1;
            shh200(:,d1) = (wtk.*hess(:,kr2,d1))*h200_1;
            shh201(:,d1) = (wtk.*hess(:,kr2,d1))*h201_1;
            shh020(:,d1) = (wtk.*hess(:,kr2,d1))*h020_1;
            shh002(:,d1) = (wtk.*hess(:,kr2,d1))*h002_1;
            shh110(:,d1) = (wtk.*hess(:,kr2,d1))*h110_1;
            shh011(:,d1) = (wtk.*hess(:,kr2,d1))*h011_1;
            shh120(:,d1) = (wtk.*hess(:,kr2,d1))*h120_1;
            shh210(:,d1) = (wtk.*hess(:,kr2,d1))*h210_1;
            shh021(:,d1) = (wtk.*hess(:,kr2,d1))*h021_1;
            for d2=lds.phases                
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
                stensv2(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v2_1;
                stensv22(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v22_1;
                for d3=lds.phases
                    stens4v1(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v1_1;
                    stens4v2(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v2_1;
                    for d4=lds.phases
                        stens5v1(:,d4,d3,d2,d1) = (wtk.*tens5(:,kr2,d1,d2,d3,d4))*v1_1;
                        stens5v2(:,d4,d3,d2,d1) = (wtk.*tens5(:,kr2,d1,d2,d3,d4))*v2_1;
                    end
                    stens5v1v1(:,d3,d2,d1) = stens5v1v1(:,d3,d2,d1)+wtk.*stens5v1(:,kr2,d3,d2,d1)*v1_1;
                    stens5v1v2(:,d3,d2,d1) = stens5v1v2(:,d3,d2,d1)+wtk.*stens5v1(:,kr2,d3,d2,d1)*v2_1;
                    stens5v2v2(:,d3,d2,d1) = stens5v2v2(:,d3,d2,d1)+wtk.*stens5v2(:,kr2,d3,d2,d1)*v2_1;
                end                
                sten4v1v1(:,d2,d1)=sten4v1v1(:,d2,d1)+wtk.*stens4v1(:,kr2,d2,d1)*v1_1;
                sten4v1v2(:,d2,d1)=sten4v1v2(:,d2,d1)+wtk.*stens4v1(:,kr2,d2,d1)*v2_1;
                sten4v2v2(:,d2,d1)=sten4v2v2(:,d2,d1)+wtk.*stens4v2(:,kr2,d2,d1)*v2_1;
                sten4v2v22(:,d2,d1)=sten4v2v22(:,d2,d1)+wtk.*stens4v2(:,kr2,d2,d1)*v22_1;
                
                stens5v1v1v1(:,d2,d1) = stens5v1v1v1(:,d2,d1)+wtk.*stens5v1v1(:,kr2,d2,d1)*v1_1;
                stens5v1v1v2(:,d2,d1) = stens5v1v1v2(:,d2,d1)+wtk.*stens5v1v1(:,kr2,d2,d1)*v2_1;
                stens5v1v2v2(:,d2,d1) = stens5v1v2v2(:,d2,d1)+wtk.*stens5v1v2(:,kr2,d2,d1)*v2_1;
                stens5v2v2v2(:,d2,d1) = stens5v2v2v2(:,d2,d1)+wtk.*stens5v2v2(:,kr2,d2,d1)*v2_1;
            end  
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
            stenv1v2(:,d1)=stenv1v2(:,d1)+wtk.*stensv1(:,kr2,d1)*v2_1;
            stenv1v22(:,d1)=stenv1v22(:,d1)+wtk.*stensv1(:,kr2,d1)*v22_1;
            stenv1h200(:,d1)=stenv1h200(:,d1)+wtk.*stensv1(:,kr2,d1)*h200_1;
            stenv1h020(:,d1)=stenv1h020(:,d1)+wtk.*stensv1(:,kr2,d1)*h020_1;
            stenv1h101(:,d1)=stenv1h101(:,d1)+wtk.*stensv1(:,kr2,d1)*h101_1;
            stenv1h011(:,d1)=stenv1h011(:,d1)+wtk.*stensv1(:,kr2,d1)*h011_1;
            stenv1h110(:,d1)=stenv1h110(:,d1)+wtk.*stensv1(:,kr2,d1)*h110_1;
            stenv2v2(:,d1)=stenv2v2(:,d1)+wtk.*stensv2(:,kr2,d1)*v2_1;
            stenv2v22(:,d1)=stenv2v22(:,d1)+wtk.*stensv2(:,kr2,d1)*v22_1;
            stenv22v22(:,d1)=stenv22v22(:,d1)+wtk.*stensv22(:,kr2,d1)*v22_1;
            stenv2h200(:,d1)=stenv2h200(:,d1)+wtk.*stensv2(:,kr2,d1)*h200_1;
            stenv2h002(:,d1)=stenv2h002(:,d1)+wtk.*stensv2(:,kr2,d1)*h002_1;
            stenv2h011(:,d1)=stenv2h011(:,d1)+wtk.*stensv2(:,kr2,d1)*h011_1;
            stenv2h101(:,d1)=stenv2h101(:,d1)+wtk.*stensv2(:,kr2,d1)*h101_1;
            stenv22h110(:,d1)=stenv22h110(:,d1)+wtk.*stensv22(:,kr2,d1)*h110_1;
            stenv22h101(:,d1)=stenv22h101(:,d1)+wtk.*stensv22(:,kr2,d1)*h101_1;
            stenv22h200(:,d1)=stenv22h200(:,d1)+wtk.*stensv22(:,kr2,d1)*h200_1;
            stenv22h020(:,d1)=stenv22h020(:,d1)+wtk.*stensv22(:,kr2,d1)*h020_1;
            
            sten4v1v1v1(:,d1)=sten4v1v1v1(:,d1)+wtk.*sten4v1v1(:,kr2,d1)*v1_1;
            sten4v1v1v2(:,d1)=sten4v1v1v2(:,d1)+wtk.*sten4v1v1(:,kr2,d1)*v2_1;
            sten4v1v1v22(:,d1)=sten4v1v1v22(:,d1)+wtk.*sten4v1v1(:,kr2,d1)*v22_1;
            sten4v1v2v2(:,d1)=sten4v1v2v2(:,d1)+wtk.*sten4v1v2(:,kr2,d1)*v2_1;
            sten4v1v2v22(:,d1)=sten4v1v2v22(:,d1)+wtk.*sten4v1v2(:,kr2,d1)*v22_1;
            sten4v2v2v2(:,d1)=sten4v2v2v2(:,d1)+wtk.*sten4v2v2(:,kr2,d1)*v2_1;
            sten4v2v2v22(:,d1)=sten4v2v2v22(:,d1)+wtk.*sten4v2v2(:,kr2,d1)*v22_1;
            sten4v2v22v22(:,d1)=sten4v2v22v22(:,d1)+wtk.*sten4v2v22(:,kr2,d1)*v22_1;
            
            stens5v1v1v1v1(:,d1) = stens5v1v1v1v1(:,d1)+wtk.*stens5v1v1v1(:,kr2,d1)*v1_1;
            stens5v1v1v1v2(:,d1) = stens5v1v1v1v2(:,d1)+wtk.*stens5v1v1v1(:,kr2,d1)*v2_1;
            stens5v1v1v2v2(:,d1) = stens5v1v1v2v2(:,d1)+wtk.*stens5v1v1v2(:,kr2,d1)*v2_1;
            stens5v1v2v2v22(:,d1) = stens5v1v2v2v22(:,d1)+wtk.*stens5v1v2v2(:,kr2,d1)*v22_1;
            stens5v2v2v2v22(:,d1) = stens5v2v2v2v22(:,d1)+wtk.*stens5v2v2v2(:,kr2,d1)*v22_1;
        end
        fxhessv1h400(range2)=wtk.*shv1(:,kr2)*h400_1;
        fxhessv1h310(range2)=wtk.*shv1(:,kr2)*h310_1;
        fxhessv1h211(range2)=wtk.*shv1(:,kr2)*h211_1;
        fxhessv1h022(range2)=wtk.*shv1(:,kr2)*h022_1;
        fxhessv1h121(range2)=wtk.*shv1(:,kr2)*h121_1;
        fxhessv2h400(range2)=wtk.*shv2(:,kr2)*h400_1;
        fxhessv2h301(range2)=wtk.*shv2(:,kr2)*h301_1;
        fxhessv2h211(range2)=wtk.*shv2(:,kr2)*h211_1;
        fxhessv2h112(range2)=wtk.*shv2(:,kr2)*h112_1;
        fxhessv2h022(range2)=wtk.*shv2(:,kr2)*h022_1;
        fxhessv22h220(range2)=wtk.*shv22(:,kr2)*h220_1;
        fxhessv22h031(range2)=wtk.*shv22(:,kr2)*h031_1;
        fxhessh200h300(range2)=wtk.*shh200(:,kr2)*h300_1;
        fxhessh200h111(range2)=wtk.*shh200(:,kr2)*h111_1;
        fxhessh200h210(range2)=wtk.*shh200(:,kr2)*h210_1;
        fxhessh200h021(range2)=wtk.*shh200(:,kr2)*h021_1;
        fxhessh201h020(range2)=wtk.*shh201(:,kr2)*h020_1;
        fxhessh020h102(range2)=wtk.*shh020(:,kr2)*h102_1;
        fxhessh020h021(range2)=wtk.*shh020(:,kr2)*h021_1;
        fxhessh002h030(range2)=wtk.*shh002(:,kr2)*h030_1;
        fxhessh110h300(range2)=wtk.*shh110(:,kr2)*h300_1;
        fxhessh110h201(range2)=wtk.*shh110(:,kr2)*h201_1;
        fxhessh110h012(range2)=wtk.*shh110(:,kr2)*h012_1;
        fxhessh110h111(range2)=wtk.*shh110(:,kr2)*h111_1;
        fxhessh011h300(range2)=wtk.*shh011(:,kr2)*h300_1;
        fxhessh011h111(range2)=wtk.*shh011(:,kr2)*h111_1;
        fxhessh120h101(range2)=wtk.*shh120(:,kr2)*h101_1;
        fxhessh210h011(range2)=wtk.*shh210(:,kr2)*h011_1;
        fxhessh021h011(range2)=wtk.*shh021(:,kr2)*h011_1;
        
        fxtensv1v1h300(range2)=wtk.*stenv1v1(:,kr2)*h300_1;
        fxtensv1v1h210(range2)=wtk.*stenv1v1(:,kr2)*h210_1;
        fxtensv1v1h021(range2)=wtk.*stenv1v1(:,kr2)*h021_1;
        fxtensv1v1h111(range2)=wtk.*stenv1v1(:,kr2)*h111_1;
        fxtensv1v2h012(range2)=wtk.*stenv1v2(:,kr2)*h012_1;
        fxtensv1v2h111(range2)=wtk.*stenv1v2(:,kr2)*h111_1;
        fxtensv1v2h300(range2)=wtk.*stenv1v2(:,kr2)*h300_1;
        fxtensv1v2h201(range2)=wtk.*stenv1v2(:,kr2)*h201_1;        
        fxtensv1v22h120(range2)=wtk.*stenv1v22(:,kr2)*h120_1;        
        fxtensv1h200h200(range2)=wtk.*stenv1h200(:,kr2)*h200_1;
        fxtensv1h200h110(range2)=wtk.*stenv1h200(:,kr2)*h110_1;
        fxtensv1h020h101(range2)=wtk.*stenv1h020(:,kr2)*h101_1;
        fxtensv1h020h002(range2)=wtk.*stenv1h020(:,kr2)*h002_1;
        fxtensv1h200h011(range2)=wtk.*stenv1h200(:,kr2)*h011_1;        
        fxtensv1h101h110(range2)=wtk.*stenv1h101(:,kr2)*h110_1;        
        fxtensv1h011h011(range2)=wtk.*stenv1h011(:,kr2)*h011_1;        
        fxtensv1h110h011(range2)=wtk.*stenv1h110(:,kr2)*h011_1;        
        fxtensv2v2h201(range2)=wtk.*stenv2v2(:,kr2)*h201_1;        
        fxtensv2v2h102(range2)=wtk.*stenv2v2(:,kr2)*h102_1;        
        fxtensv2v2h012(range2)=wtk.*stenv2v2(:,kr2)*h012_1;        
        fxtensv2v22h300(range2)=wtk.*stenv2v22(:,kr2)*h300_1;        
        fxtensv2v22h210(range2)=wtk.*stenv2v22(:,kr2)*h210_1;        
        fxtensv2v22h111(range2)=wtk.*stenv2v22(:,kr2)*h111_1;        
        fxtensv2v22h021(range2)=wtk.*stenv2v22(:,kr2)*h021_1;                
        fxtensv22v22h030(range2)=wtk.*stenv22v22(:,kr2)*h030_1;                
        fxtensv2h200h200(range2)=wtk.*stenv2h200(:,kr2)*h200_1;
        fxtensv2h002h020(range2)=wtk.*stenv2h002(:,kr2)*h020_1;
        fxtensv2h200h101(range2)=wtk.*stenv2h200(:,kr2)*h101_1;
        fxtensv2h101h110(range2)=wtk.*stenv2h101(:,kr2)*h110_1;
        fxtensv2h200h011(range2)=wtk.*stenv2h200(:,kr2)*h011_1;
        fxtensv2h002h110(range2)=wtk.*stenv2h002(:,kr2)*h110_1;
        fxtensv2h011h101(range2)=wtk.*stenv2h011(:,kr2)*h101_1;
        fxtensv2h011h011(range2)=wtk.*stenv2h011(:,kr2)*h011_1;
        fxtensv22h110h110(range2)=wtk.*stenv22h110(:,kr2)*h110_1;        
        fxtensv22h101h110(range2)=wtk.*stenv22h101(:,kr2)*h110_1;        
        fxtensv22h200h020(range2)=wtk.*stenv22h200(:,kr2)*h020_1;        
        fxtensv22h020h011(range2)=wtk.*stenv22h020(:,kr2)*h011_1;        
        
        fxtens4v1v1v1h200(range2)=wtk.*sten4v1v1v1(:,kr2)*h200_1;
        fxtens4v1v1v1h110(range2)=wtk.*sten4v1v1v1(:,kr2)*h110_1;
        fxtens4v1v1v1h011(range2)=wtk.*sten4v1v1v1(:,kr2)*h011_1;
        fxtens4v1v1v2h200(range2)=wtk.*sten4v1v1v2(:,kr2)*h200_1;
        fxtens4v1v1v2h101(range2)=wtk.*sten4v1v1v2(:,kr2)*h101_1;
        fxtens4v1v1v2h011(range2)=wtk.*sten4v1v1v2(:,kr2)*h011_1;
        fxtens4v1v1v22h020(range2)=wtk.*sten4v1v1v22(:,kr2)*h020_1;
        fxtens4v1v2v2h101(range2)=wtk.*sten4v1v2v2(:,kr2)*h101_1;
        fxtens4v1v2v2h002(range2)=wtk.*sten4v1v2v2(:,kr2)*h002_1;
        fxtens4v1v2v22h200(range2)=wtk.*sten4v1v2v22(:,kr2)*h200_1;
        fxtens4v1v2v22h110(range2)=wtk.*sten4v1v2v22(:,kr2)*h110_1;
        fxtens4v1v2v22h011(range2)=wtk.*sten4v1v2v22(:,kr2)*h011_1;
        fxtens4v2v2v2h002(range2)=wtk.*sten4v2v2v2(:,kr2)*h002_1;
        fxtens4v2v2v22h200(range2)=wtk.*sten4v2v2v22(:,kr2)*h200_1;
        fxtens4v2v2v22h101(range2)=wtk.*sten4v2v2v22(:,kr2)*h101_1;
        fxtens4v2v2v22h011(range2)=wtk.*sten4v2v2v22(:,kr2)*h011_1;
        fxtens4v2v22v22h020(range2)=wtk.*sten4v2v22v22(:,kr2)*h020_1;
        
        fxtens5v1v1v1v1v1(range2)=wtk.*stens5v1v1v1v1(:,kr2)*v1_1;        
        fxtens5v1v1v1v1v2(range2)=wtk.*stens5v1v1v1v1(:,kr2)*v2_1;        
        fxtens5v1v1v1v2v22(range2)=wtk.*stens5v1v1v1v2(:,kr2)*v22_1;        
        fxtens5v1v1v2v2v22(range2)=wtk.*stens5v1v1v2v2(:,kr2)*v22_1;        
        fxtens5v1v2v2v22v22(range2)=wtk.*stens5v1v2v2v22(:,kr2)*v22_1;        
        fxtens5v2v2v2v22v22(range2)=wtk.*stens5v2v2v2v22(:,kr2)*v22_1;        
        range2 = range2+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
end

%computation of a500
diffh300 = T*fxjac*h300+diffh300.';
a500 = 1/(120*T^2)*w1*(fxtens5v1v1v1v1v1.'+10*fxtens4v1v1v1h200.'+10*fxtensv1v1h300.'+15*fxtensv1h200h200.'...
    +10*fxhessh200h300.'+5*fxhessv1h400.'-20*alpha2001*diffh300-120*alpha4001*fxjac*v1);

%computation of b410
diffh210 = T*fxjac*h210-i*theta*h210C+diffh210.';
b410 = 1/(24*T^2)*w2*(fxtens5v1v1v1v1v2.'+6*fxtens4v1v1v2h200.'+4*fxtens4v1v1v1h110.'+4*fxtensv1v2h300.'...
    +6*fxtensv1v1h210.'+3*fxtensv2h200h200.'+12*fxtensv1h200h110.'+4*fxhessv1h310.'+4*fxhessh110h300.'...
    +6*fxhessh200h210.'+fxhessv2h400.'-24*alpha4001*fxjac*v2-12*alpha2001*diffh210)+1/T^2*alpha4001*i*theta/T;

%computation of a311
diffh111 = T*fxjac*h111+diffh111.';
a311 = 1/(6*T^2)*w1*(real(fxtens5v1v1v1v2v22)'+fxtens4v1v1v1h011'+6*real(fxtens4v1v1v2h101)'+3*real(fxtens4v1v2v22h200)'...
    +3*fxtensv1h200h011'+6*real(fxtensv2h200h101)'+6*real(fxtensv1v2h201)'+real(fxtensv2v22h300)'+3*fxtensv1v1h111'...
    +6*real(fxtensv1h101h110)'+3*fxhessh200h111'+6*real(fxhessh110h201)'+2*real(fxhessv2h301)'+fxhessh011h300'...
    +3*fxhessv1h211'-6*alpha2111*fxjac*v1-alpha0111*diffh300-6*alpha2001*diffh111);

%computation of b221
diffh021 = T*fxjac*h021-i*theta*h021C+diffh021.';
b221 = 1/(4*T^2)*w2*(fxtens5v1v1v2v2v22.'+fxtens4v2v2v22h200.'+2*fxtens4v1v2v2h101.'+2*fxtens4v1v1v2h011.'...
    +fxtens4v1v1v22h020.'+fxtens4v1v2v22h110.'+2*fxtensv22h110h110.'+fxtensv1v1h021.'+fxtensv2v2h201.'...
    +fxtensv22h200h020.'+2*fxtensv2v22h210.'+2*fxtensv1v22h120.'+2*fxtensv1h020h101.'+4*fxtensv1v2h111.'...
    +4*fxtensv2h101h110.'+2*fxtensv2h200h011.'+4*fxtensv1h110h011.'+fxhessv22h220.'+2*fxhessv1h121.'...
    +2*fxhessh120h101.'+4*fxhessh110h111.'+2*fxhessh210h011.'+2*fxhessv2h211.'+fxhessh200h021.'+fxhessh201h020.'...
    -2*alpha0111*diffh210-4*alpha2111*fxjac*v2-2*alpha2001*diffh021)+i*theta/T^3*alpha2111;

%computation of a122
a122 = 1/(4*T^2)*w1*(real(fxtens5v1v2v2v22v22)'+4*real(fxtens4v2v2v22h101)'+2*real(fxtens4v1v2v2h002)'...
    +4*real(fxtens4v1v2v22h011)'+2*fxtensv1h011h011'+8*real(fxtensv2h011h101)'+4*real(fxtensv2h002h110)'...
    +real(fxtensv1h020h002)'+4*real(fxtensv1v2h012)'+4*real(fxtensv2v22h111)'+2*real(fxtensv2v2h102)'...
    +4*fxhessh011h111'+2*real(fxhessh020h102)'+4*real(fxhessv2h112)'+fxhessv1h022'+4*real(fxhessh110h012)'...
    -4*alpha0111*diffh111-4*alpha0221*fxjac*v1);

%computation of b032
b032 = 1/(12*T^2)*w2*(fxtens5v2v2v2v22v22.'+fxtens4v2v2v2h002.'+3*fxtens4v2v22v22h020.'+6*fxtens4v2v2v22h011.'...
    +6*fxtensv22h020h011.'+6*fxtensv2v22h021.'+fxtensv22v22h030.'+3*fxtensv2h002h020.'+3*fxtensv2v2h012.'...
    +6*fxtensv2h011h011.'+3*fxhessh020h021.'+6*fxhessh021h011.'+3*fxhessv2h022.'+2*fxhessv22h031.'+fxhessh002h030.'...
    -6*alpha0111*diffh021-12*alpha0221*fxjac*v2)+alpha0221*i*theta/T^3;
end

%transformation formulas
G2100 = a300;
G1011 = a111;
H1110 = b210-1i*theta/T*alpha200;
H0021 = b021-1i*theta/T/alpha011;

if higherorderderivatives == 1
G3200 = a500-a300*alpha200;
G2111 = -a111*alpha200-a300*alpha011+a311;
G1022 = a122-a111*alpha011;
H2210 = b410-1i*theta/T*alpha400+(-b210+1i*theta/T*alpha200)*alpha200;
H1121 = -1i*theta/T*alpha211+b221+(-b021+1i*theta/T*alpha011)*alpha200+(-b210+1i*theta/T*alpha200)*alpha011;
H0032 = b032-1i*theta/T*alpha022+(-b021+1i*theta/T*alpha011)*alpha011;
end 

p11 = real(G2100);
p12 = real(G1011);
p21 = real(H1110);
p22 = real(H0021);
theta = p12/p22;
delta = p21/p11;
if higherorderderivatives == 1
s1 = real(G1022)+real(G1011)*( real(H1121)/real(H1110)-2*real(H0032)/real(H0021)-real(G3200)*real(H0021)/(real(G2100)*real(H1110)) );
s2 = real(H2210)+real(H1110)*( real(G2111)/real(G1011)-2*real(G3200)/real(G2100)-real(G2100)*real(H0032)/(real(G1011)*real(H0021)) );
theta_big = s1/p22^2;
delta_big = s2/p11^2;
l1 = (delta-1)/(theta-1)*( theta*(theta-1)*delta_big+delta*(delta-1)*theta_big );
else
l1 = NaN;
end

result = [p11 p22 theta delta sign(l1)];

% ---------------------------------------------------------------
function [x,p,T] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global lds

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+1:lds.PeriodIdx+2);
x = x0(lds.coords);
T = x0(lds.PeriodIdx);

% -------------------------------------------------------------

function [jac1,jac2] = BVP_jac(BVP_jac_f,BVP_jac_bc1,x,p,T,theta,pars,nc)
global lds 

ups = reshape(x,lds.nphase,lds.tps);
p = num2cell(p);

jac1 = spalloc(lds.ncoords,lds.ncoords,(lds.ncol+4)*lds.nphase);

range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
jac2 = jac1;
eitheta = 1i*theta*lds.pwwt;
for jk = lds.tsts
  % value of polynomial on each collocation point
  xp = ups(:,range0)*lds.wt;

  % evaluate part of Jacobian matrix
  jac2(range1,range2) = feval(BVP_jac_f,lds.func,xp,p,T,jk);
  jac1(range1,range2) = jac2(range1,range2)+eitheta;

  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end

% boundary conditions
range = (lds.tps-1)*lds.nphase + (lds.phases);
jac1(range,[lds.phases range]) = feval(BVP_jac_bc1);
jac2(range,[lds.phases range]) = jac1(range,[lds.phases range]);