function result=nf_NSNS(x)

% calculates NSNS normal form coefficient
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
                idx1 = jk;
                idx2 = jk+idx;
                smallest_sum = val;
            end
        end
    end
theta1 = abs(angle(d(idx2)));

%idsx1 is smaller than idx2
d(idx2) = [];
d(idx1) = [];

smallest_sum = Inf;
    for jk=1:lds.nphase-3
        [val,idx] = min(abs(d(jk+1:lds.nphase-2)*d(jk)-1));
        if val < smallest_sum
            idx1 = jk;
            idx2 = jk+idx;
            smallest_sum = val;
        end
    end
theta2 = abs(angle(d(idx2)));


%computation main part matrices
[Jjac1,Jjac2,Jjac3] = BVP_jac('bordBVP_R2_f','BVP_LC1_jac_bc',xt,p,T,theta1,theta2,1,1);

%computation of v1
b = []; b(lds.ncoords+1)=1; b=b';
jace=[Jjac1,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

Jjac1 = [Jjac1 p;q' 0];
[LJ,UJ] = lu(Jjac1);
vext = UJ \(LJ\b);

%rescale vext
v1 = vext(1:lds.ncoords);
vps = reshape(v1,lds.nphase,lds.tps);
ficd = dot(vps,vps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v1 = v1/sqrt(ic);

%computation of v1*
wext = LJ'\(UJ'\b);
w1 = wext(1:lds.ncoords-lds.nphase);

%computation of v2
b = []; b(lds.ncoords+1)=1; b=b';
jace=[Jjac2,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

Jjac2 = [Jjac2 p;q' 0];
[LJ,UJ] = lu(Jjac2);
vext = UJ \(LJ\b);

%rescale vext
v2 = vext(1:lds.ncoords);
vps = reshape(v2,lds.nphase,lds.tps);
ficd = dot(vps,vps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v2 = v2/sqrt(ic);

%computation of v2*
wext = LJ'\(UJ'\b);
w2 = wext(1:lds.ncoords-lds.nphase);

%computation of phi*
jace=[Jjac3,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

Jjac3=[Jjac3 p;q' 0];

%computation of phi*
phi = Jjac3'\b;
phi = phi(1:lds.ncoords-lds.nphase);

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
phi = phi/(phi'*f1);

%rescale v1*
w1 = w1/(w1'*v1C)';

%rescale v2*
w2 = w2/(w2'*v2C)';

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
w11ps = conj(w1ps);
w22ps = conj(w2ps);
icphi = zeros(1,lds.ncoords);
icw11 = zeros(1,lds.ncoords);
icw22 = zeros(1,lds.ncoords);
wt = lds.wt';
for j=lds.tsts
    % value of polynomial on each collocation point
    xpp(:,range6) = ups(:,range1)*lds.wt;
    xp = xpp(:,range6); 
    v1_1 = v1(range4);
    v11_1 = conj(v1_1);
    v2_1 = v2(range4);
    v22_1 = conj(v2_1);
    p = phips(:,range6)*wt;
    icphi(range4) = icphi(range4)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    p = w11ps(:,range6)*wt;
    icw11(range4) = icw11(range4)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    p = w22ps(:,range6)*wt;
    icw22(range4) = icw22(range4)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
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
        fxhessv1v11(range2)=wtk.*shv1(:,kr2)*v11_1;
        fxhessv1v2(range2)=wtk.*shv1(:,kr2)*v2_1;
        fxhessv1v22(range2)=wtk.*shv1(:,kr2)*v22_1;
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

range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
range3 = lds.cols;
twitheta1 = 2*1i*theta1*lds.pwwt;
twitheta2 = 2*1i*theta2*lds.pwwt;
theta1plustheta2plus = 1i*theta1*lds.pwwt+1i*theta2*lds.pwwt;
theta1plustheta2min = 1i*theta1*lds.pwwt-1i*theta2*lds.pwwt;
theta1mintheta2plus = -1i*theta1*lds.pwwt+1i*theta2*lds.pwwt;
threeitheta1 = 3*1i*theta1*lds.pwwt;
threeitheta2 = 3*1i*theta2*lds.pwwt;
twotheta1plustheta2plus = 2i*theta1*lds.pwwt+1i*theta2*lds.pwwt;
twotheta1plustheta2min = 2i*theta1*lds.pwwt-1i*theta2*lds.pwwt;
theta1min2theta2plus = -1i*theta1*lds.pwwt+2i*theta2*lds.pwwt;
theta1plus2theta2plus = 1i*theta1*lds.pwwt+2i*theta2*lds.pwwt;
fourtheta1plus = 4i*theta1*lds.pwwt;
threetheta1plustheta2plus = 3*i*theta1*lds.pwwt+i*theta2*lds.pwwt;
threetheta1plustheta2min = 3*i*theta1*lds.pwwt-i*theta2*lds.pwwt;
twotheta1plus2theta2plus = 2*i*theta1*lds.pwwt+2*i*theta2*lds.pwwt;
twotheta1plus2theta2min = 2*i*theta1*lds.pwwt-2*i*theta2*lds.pwwt;

J2theta1plus = Jjac1;
J2theta2plus = Jjac1;
Jtheta1plustheta2plus = Jjac1;
Jtheta1plustheta2min = Jjac1;
Jtheta1mintheta2plus = Jjac1;
J3theta1plus = Jjac1;
J3theta2plus = Jjac1;
J2theta1plustheta2plus = Jjac1;
J2theta1plustheta2min = Jjac1;
Jtheta1min2theta2plus = Jjac1;
Jtheta1plus2theta2plus = Jjac1;
J4theta1plus = Jjac1;
J3theta1plustheta2plus = Jjac1;
J3theta1plustheta2min = Jjac1;
J2theta1plus2theta2plus = Jjac1;
J2theta1plus2theta2min = Jjac1;

for jk  = lds.tsts
  % evaluate part of Jacobian matrix
  J2theta1plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+twitheta1;  
  J2theta2plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+twitheta2;  
  Jtheta1plustheta2plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+theta1plustheta2plus;  
  Jtheta1plustheta2min(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+theta1plustheta2min;  
  Jtheta1mintheta2plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+theta1mintheta2plus;  
  J3theta1plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+threeitheta1;  
  J3theta2plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+threeitheta2;  
  J2theta1plustheta2plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+twotheta1plustheta2plus;  
  J2theta1plustheta2min(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+twotheta1plustheta2min;  
  Jtheta1min2theta2plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+theta1min2theta2plus;  
  Jtheta1plus2theta2plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+theta1plus2theta2plus;  
  J4theta1plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+fourtheta1plus;  
  J3theta1plustheta2plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+threetheta1plustheta2plus;  
  J3theta1plustheta2min(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+threetheta1plustheta2min;
  J2theta1plus2theta2plus(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+twotheta1plus2theta2plus;  
  J2theta1plus2theta2min(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+twotheta1plus2theta2min;  
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
  range3 = range3 + lds.ncol;
end
% remove borders
J2theta1plus(:,end)=[];J2theta1plus(end,:)=[];
J2theta2plus(:,end)=[];J2theta2plus(end,:)=[];
Jtheta1plustheta2plus(:,end)=[];Jtheta1plustheta2plus(end,:)=[];
Jtheta1plustheta2min(:,end)=[];Jtheta1plustheta2min(end,:)=[];
Jtheta1mintheta2plus(:,end)=[];Jtheta1mintheta2plus(end,:)=[];
J3theta1plus(:,end)=[];J3theta1plus(end,:)=[];
J3theta2plus(:,end)=[];J3theta2plus(end,:)=[];
J2theta1plustheta2plus(:,end)=[];J2theta1plustheta2plus(end,:)=[];
J2theta1plustheta2min(:,end)=[];J2theta1plustheta2min(end,:)=[];
Jtheta1min2theta2plus(:,end)=[];Jtheta1min2theta2plus(end,:)=[];
Jtheta1plus2theta2plus(:,end)=[];Jtheta1plus2theta2plus(end,:)=[];
J4theta1plus(:,end)=[];J4theta1plus(end,:)=[];
J3theta1plustheta2plus(:,end)=[];J3theta1plustheta2plus(end,:)=[];
J3theta1plustheta2min(:,end)=[];J3theta1plustheta2min(end,:)=[];
J2theta1plus2theta2plus(:,end)=[];J2theta1plus2theta2plus(end,:)=[];
J2theta1plus2theta2min(:,end)=[];J2theta1plus2theta2min(end,:)=[];

%computation of h2000 and h0020
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v1.'; rl(lds.ncoords) = 0;
diffh2000 = rl(1:lds.ncoords-lds.nphase);
h2000 = J2theta1plus\rl.';
rl(1:lds.ncoords-lds.nphase) = T*fxhessv2v2.'; rl(lds.ncoords) = 0;
diffh0020 = rl(1:lds.ncoords-lds.nphase);
h0020 = J2theta2plus\rl.';

%computation of alpha11001 and alpha00111
alpha11001 = phi'*fxhessv1v11.';
alpha00111 = phi'*fxhessv2v22.';
alpha1100 = alpha11001/T;
alpha0011 = alpha00111/T;

%computation of h1100 and h0011
Jjac3(lds.ncoords+1,lds.coords) = icphi;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxhessv1v11)'-T*alpha11001*f1; rl(lds.ncoords+1) = 0;
diffh1100 = rl(1:lds.ncoords-lds.nphase);
h1100 = Jjac3\rl.';
h1100 = h1100(lds.coords);

rl(1:lds.ncoords-lds.nphase) = T*real(fxhessv2v22)'-T*alpha00111*f1; rl(lds.ncoords+1) = 0;
diffh0011 = rl(1:lds.ncoords-lds.nphase);
h0011 = Jjac3\rl.';
h0011 = h0011(lds.coords);

%computation of h1010
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v2.'; rl(lds.ncoords) = 0;
diffh1010 = rl(1:lds.ncoords-lds.nphase);
h1010 = Jtheta1plustheta2plus\rl.';

%computation of h1001
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v22.'; rl(lds.ncoords) = 0;
diffh1001 = rl(1:lds.ncoords-lds.nphase);
h1001 = Jtheta1plustheta2min\rl.';

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
    v1_1 = v1(range4);
    v11_1 = conj(v1_1);
    v2_1 = v2(range4);
    v22_1 = conj(v2_1);
    h1100_1 = h1100(range4);
    h0011_1 = h0011(range4);
    h1010_1 = h1010(range4);
    h1001_1 = h1001(range4);
    h2000_1 = h2000(range4);
    h0110_1 = conj(h1001_1);
    h0020_1 = h0020(range4);
    h0101_1 = conj(h1010_1);
    range5 = lds.phases;         
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        stenv1v1 = zeros(lds.nphase,lds.nphase);
        stenv1v11 = zeros(lds.nphase,lds.nphase);
        stenv1v2 = zeros(lds.nphase,lds.nphase);
        stenv11v2 = zeros(lds.nphase,lds.nphase);
        stenv2v2 = zeros(lds.nphase,lds.nphase);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shv11(:,d1) = (wtk.*hess(:,kr2,d1))*v11_1;
            shv2(:,d1) = (wtk.*hess(:,kr2,d1))*v2_1;
            shv22(:,d1) = (wtk.*hess(:,kr2,d1))*v22_1;
            for d2=lds.phases
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
                stensv11(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v11_1;
                stensv2(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v2_1;
            end    
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
            stenv1v11(:,d1)=stenv1v11(:,d1)+wtk.*stensv1(:,kr2,d1)*v11_1;
            stenv1v2(:,d1)=stenv1v2(:,d1)+wtk.*stensv1(:,kr2,d1)*v2_1;
            stenv11v2(:,d1)=stenv11v2(:,d1)+wtk.*stensv11(:,kr2,d1)*v2_1;
            stenv2v2(:,d1)=stenv2v2(:,d1)+wtk.*stensv2(:,kr2,d1)*v2_1;
        end
        fxhessv1h2000(range2)=wtk.*shv1(:,kr2)*h2000_1;
        fxhessv1h1100(range2)=wtk.*shv1(:,kr2)*h1100_1;
        fxhessv1h1001(range2)=wtk.*shv1(:,kr2)*h1001_1;
        fxhessv1h1010(range2)=wtk.*shv1(:,kr2)*h1010_1;
        fxhessv1h0110(range2)=wtk.*shv1(:,kr2)*h0110_1;
        fxhessv1h0011(range2)=wtk.*shv1(:,kr2)*h0011_1;
        fxhessv1h0101(range2)=wtk.*shv1(:,kr2)*h0101_1;
        fxhessv1h0020(range2)=wtk.*shv1(:,kr2)*h0020_1;
        fxhessv11h2000(range2)=wtk.*shv11(:,kr2)*h2000_1;
        fxhessv11h1001(range2)=wtk.*shv11(:,kr2)*h1001_1;
        fxhessv11h0020(range2)=wtk.*shv11(:,kr2)*h0020_1;
        fxhessv11h1010(range2)=wtk.*shv11(:,kr2)*h1010_1;
        fxhessv11h0011(range2)=wtk.*shv11(:,kr2)*h0011_1;
        fxhessv2h2000(range2)=wtk.*shv2(:,kr2)*h2000_1;
        fxhessv2h0110(range2)=wtk.*shv2(:,kr2)*h0110_1;
        fxhessv2h0011(range2)=wtk.*shv2(:,kr2)*h0011_1;
        fxhessv2h1100(range2)=wtk.*shv2(:,kr2)*h1100_1;
        fxhessv2h1010(range2)=wtk.*shv2(:,kr2)*h1010_1;
        fxhessv2h1001(range2)=wtk.*shv2(:,kr2)*h1001_1;
        fxhessv2h0020(range2)=wtk.*shv2(:,kr2)*h0020_1;        
        fxhessv2h0101(range2)=wtk.*shv2(:,kr2)*h0101_1;
        fxhessv22h0110(range2)=wtk.*shv22(:,kr2)*h0110_1;
        fxhessv22h2000(range2)=wtk.*shv22(:,kr2)*h2000_1;
        fxhessv22h0020(range2)=wtk.*shv22(:,kr2)*h0020_1;
        fxhessv22h1010(range2)=wtk.*shv22(:,kr2)*h1010_1;
        fxhessv22h1100(range2)=wtk.*shv22(:,kr2)*h1100_1;
        
        fxtensv1v1v1(range2)=wtk.*stenv1v1(:,kr2)*v1_1;
        fxtensv1v1v11(range2)=wtk.*stenv1v1(:,kr2)*v11_1;
        fxtensv1v1v2(range2)=wtk.*stenv1v1(:,kr2)*v2_1;
        fxtensv1v1v22(range2)=wtk.*stenv1v1(:,kr2)*v22_1;
        fxtensv1v11v2(range2)=wtk.*stenv1v11(:,kr2)*v2_1;
        fxtensv1v11v22(range2)=wtk.*stenv1v11(:,kr2)*v22_1;
        fxtensv1v2v2(range2)=wtk.*stenv1v2(:,kr2)*v2_1;
        fxtensv1v2v22(range2)=wtk.*stenv1v2(:,kr2)*v22_1;
        fxtensv11v2v2(range2)=wtk.*stenv11v2(:,kr2)*v2_1;
        fxtensv11v2v22(range2)=wtk.*stenv11v2(:,kr2)*v22_1;
        fxtensv2v2v2(range2)=wtk.*stenv2v2(:,kr2)*v2_1;
        fxtensv2v2v22(range2)=wtk.*stenv2v2(:,kr2)*v22_1;

        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
end

%computation of h3000
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v1v1.'+3*T*fxhessv1h2000.'; rl(lds.ncoords) = 0;
h3000 = J3theta1plus\rl.';

%computation of h0030
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv2v2v2.'+3*T*fxhessv2h0020.'; rl(lds.ncoords) = 0;
h0030 = J3theta2plus\rl.';

%computation of a2100
a21001 = 1/2*w1'*(fxtensv1v1v11.'+2*fxhessv1h1100.'+fxhessv11h2000.'-2*alpha11001*fxjac*v1)+i*theta1/T*alpha11001;
a2100 = a21001/T;

%computation of h2100
Jjac1(lds.ncoords+1,lds.coords) = icw11;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v1v11.'+2*T*fxhessv1h1100.'+T*fxhessv11h2000.'-T*2*a21001*v1C...
    -2*T*alpha11001*fxjac*v1+2*T*alpha11001*i*theta1/T*v1C; rl(lds.ncoords+1) = 0;
diffh2100 = rl(1:lds.ncoords-lds.nphase);
h2100 = Jjac1\rl.';
h2100 = h2100(lds.coords);

%computation of h2010
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v1v2.'+T*fxhessv2h2000.'+2*T*fxhessv1h1010.'; rl(lds.ncoords) = 0;
h2010 = J2theta1plustheta2plus\rl.';

%computation of h2001
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v1v22.'+T*fxhessv22h2000.'+2*T*fxhessv1h1001.'; rl(lds.ncoords) = 0;
h2001 = J2theta1plustheta2min\rl.';

%computation of h1020
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v2v2.'+T*fxhessv1h0020.'+2*T*fxhessv2h1010.'; rl(lds.ncoords) = 0;
h1020 = Jtheta1plus2theta2plus\rl.';

%computation of h0120
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv11v2v2.'+T*fxhessv11h0020.'+2*T*fxhessv2h0110.'; rl(lds.ncoords) = 0;
h0120 = Jtheta1min2theta2plus\rl.';

%computation of b0021
b00211 = 1/2*w2'*(fxtensv2v2v22.'+fxhessv22h0020.'+2*fxhessv2h0011.'-2*alpha00111*fxjac*v2)+i*theta2/T*alpha00111;
b0021 = b00211/T;

%computation of h0021
Jjac2(lds.ncoords+1,lds.coords) = icw22;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv2v2v22.'+T*fxhessv22h0020.'+2*T*fxhessv2h0011.'-T*2*b00211*v2C...
    -2*T*alpha00111*fxjac*v2+2*T*alpha00111*i*theta2/T*v2C; rl(lds.ncoords+1) = 0;
diffh0021 = rl(1:lds.ncoords-lds.nphase);
h0021 = Jjac2\rl.';
h0021 = h0021(lds.coords);

%computation of b1110
b11101 = w2'*(fxtensv1v11v2.'+fxhessv1h0110.'+fxhessv11h1010.'+fxhessv2h1100.'-alpha11001*fxjac*v2)+i*theta2/T*alpha11001;
b1110 = b11101/T;

%computation of b11011
b11011 = w2.'*(fxtensv1v11v22.'+fxhessv1h0101.'+fxhessv11h1001.'+fxhessv22h1100.'-alpha11001*fxjac*conj(v2))-i*theta2/T*alpha11001;
b1101 = b11011/T;

%computation of b01111
b01111 = w1.'*(fxtensv11v2v22.'+fxhessv2h0101.'+fxhessv22h0110.'+fxhessv11h0011.'-alpha00111*fxjac*conj(v1))-i*theta1/T*alpha00111;
b0111 = b01111/T;

%computation of h1110
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v11v2.'+T*fxhessv1h0110.'+T*fxhessv11h1010.'+T*fxhessv2h1100.'...
    -T*b11101*v2C-T*alpha11001*fxjac*v2+T*alpha11001*i*theta2/T*v2C; rl(lds.ncoords+1) = 0;
diffh1110 = rl(1:lds.ncoords-lds.nphase);
h1110 = Jjac2\rl.';
h1110 = h1110(lds.coords);

%computation of a10111
a10111 = w1'*(fxtensv1v2v22.'+fxhessv1h0011.'+fxhessv2h1001.'+fxhessv22h1010.'-alpha00111*fxjac*v1)+i*theta1/T*alpha00111;
a1011 = a10111/T;

%computation of a01111
a01111 = conj(w1)'*(fxtensv11v2v22.'+fxhessv11h0011.'+fxhessv2h0101.'+fxhessv22h0110.'-alpha00111*fxjac*conj(v1))-i*theta1/T*alpha00111;
a0111 = a01111/T;

%computation of h1011
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v2v22.'+T*fxhessv1h0011.'+T*fxhessv2h1001.'+T*fxhessv22h1010.'...
    -T*a10111*v1C-T*alpha00111*fxjac*v1+T*alpha00111*i*theta1/T*v1C; rl(lds.ncoords+1) = 0;
diffh1011 = rl(1:lds.ncoords-lds.nphase);
h1011 = Jjac1\rl.';
h1011 = h1011(lds.coords);

%computation of h1011
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv11v2v22.'+T*fxhessv11h0011.'+T*fxhessv2h0101.'+T*fxhessv22h0110.'...
    -T*a01111*conj(v1C)-T*alpha00111*fxjac*conj(v1)-T*alpha00111*i*theta1/T*conj(v1C); rl(lds.ncoords+1) = 0;
h0111 = conj(Jjac1)\rl.';
h0111 = h0111(lds.coords);

%computation of h1011
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v11v22.'+T*fxhessv22h1100.'+T*fxhessv1h0101.'+T*fxhessv11h1001.'...
    -T*b11011*conj(v2C)-T*alpha11001*fxjac*conj(v2)-T*alpha11001*i*theta2/T*conj(v2C); rl(lds.ncoords+1) = 0;
h1101 = conj(Jjac2)\rl.';
h1101 = h1101(lds.coords);

if higherorderderivatives == 1
% function
range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
h2000ps = reshape(h2000,lds.nphase,lds.tps);
h0020ps = reshape(h0020,lds.nphase,lds.tps);
h1100ps = reshape(h1100,lds.nphase,lds.tps);
h1010ps = reshape(h1010,lds.nphase,lds.tps);
h1011ps = reshape(h1011,lds.nphase,lds.tps);
h1001ps = reshape(h1001,lds.nphase,lds.tps);
h0011ps = reshape(h0011,lds.nphase,lds.tps);
h2100ps = reshape(h2100,lds.nphase,lds.tps);
h0021ps = reshape(h0021,lds.nphase,lds.tps);
h1110ps = reshape(h1110,lds.nphase,lds.tps);
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    h2000Cps(:,range6) = h2000ps(:,range1)*lds.wt;
    h0020Cps(:,range6) = h0020ps(:,range1)*lds.wt;
    h1100Cps(:,range6) = h1100ps(:,range1)*lds.wt;
    h1010Cps(:,range6) = h1010ps(:,range1)*lds.wt;
    h1011Cps(:,range6) = h1011ps(:,range1)*lds.wt;
    h1001Cps(:,range6) = h1001ps(:,range1)*lds.wt;
    h0011Cps(:,range6) = h0011ps(:,range1)*lds.wt;
    h2100Cps(:,range6) = h2100ps(:,range1)*lds.wt;
    h0021Cps(:,range6) = h0021ps(:,range1)*lds.wt;
    h1110Cps(:,range6) = h1110ps(:,range1)*lds.wt;
        
    v1_1 = v1(range4);
    v11_1 = conj(v1_1);
    v2_1 = v2(range4);
    v22_1 = conj(v2_1);
    h0011_1 = h0011(range4);
    h0020_1 = h0020(range4);
    h0002_1 = conj(h0020_1);
    h0030_1 = h0030(range4);
    h1001_1 = h1001(range4);
    h0110_1 = conj(h1001_1);
    h1010_1 = h1010(range4);    
    h0101_1 = conj(h1010_1);    
    h1011_1 = h1011(range4);    
    h0111_1 = h0111(range4);    
    h2000_1 = h2000(range4);
    h0200_1 = conj(h2000_1);    
    h1020_1 = h1020(range4);        
    h0120_1 = h0120(range4);
    h1002_1 = conj(h0120_1);       
    h1100_1 = h1100(range4);    
    h1110_1 = h1110(range4);    
    h1101_1 = h1101(range4);        
    h2100_1 = h2100(range4);    
    h0021_1 = h0021(range4);    
    h0012_1 = conj(h0021_1);
    h1200_1 = conj(h2100_1);
    h2010_1 = h2010(range4);     
    h2001_1 = h2001(range4);      
    h3000_1 = h3000(range4);
    
    range5 = lds.phases;         
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        stens4v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v1v11 = zeros(lds.nphase,lds.nphase,lds.nphase);        
        stens4v2v2 = zeros(lds.nphase,lds.nphase,lds.nphase);  
        stens4v1v2 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v11v2 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stenv1v1 = zeros(lds.nphase,lds.nphase);
        stenv1v11 = zeros(lds.nphase,lds.nphase);
        stenv1v2 = zeros(lds.nphase,lds.nphase);
        stenv1v22 = zeros(lds.nphase,lds.nphase);
        stenv11v22 = zeros(lds.nphase,lds.nphase);
        stenv11v2 = zeros(lds.nphase,lds.nphase);
        stenv2v2 = zeros(lds.nphase,lds.nphase);
        stenv2v22 = zeros(lds.nphase,lds.nphase);
        stenv22v22 = zeros(lds.nphase,lds.nphase);
        stens4v1v1v1 = zeros(lds.nphase,lds.nphase);
        stens4v1v1v11 = zeros(lds.nphase,lds.nphase);
        stens4v1v1v2 = zeros(lds.nphase,lds.nphase);
        stens4v1v11v2 = zeros(lds.nphase,lds.nphase);
        stens4v1v1v22 = zeros(lds.nphase,lds.nphase);
        stens4v2v2v2 = zeros(lds.nphase,lds.nphase);
        stens4v2v2v22 = zeros(lds.nphase,lds.nphase);
        stens4v1v2v2 = zeros(lds.nphase,lds.nphase);
        stens4v11v2v2 = zeros(lds.nphase,lds.nphase);
        
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        tens4 = ctens4(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,xt,pt,lds.ActiveParams);
        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shv11(:,d1) = (wtk.*hess(:,kr2,d1))*v11_1;
            shv2(:,d1) = (wtk.*hess(:,kr2,d1))*v2_1;
            shv22(:,d1) = (wtk.*hess(:,kr2,d1))*v22_1;
            shh2000(:,d1) = (wtk.*hess(:,kr2,d1))*h2000_1;
            shh1100(:,d1) = (wtk.*hess(:,kr2,d1))*h1100_1;
            shh1010(:,d1) = (wtk.*hess(:,kr2,d1))*h1010_1;
            shh1001(:,d1) = (wtk.*hess(:,kr2,d1))*h1001_1;
            shh0110(:,d1) = (wtk.*hess(:,kr2,d1))*h0110_1;
            shh0101(:,d1) = (wtk.*hess(:,kr2,d1))*h0101_1;
            shh0011(:,d1) = (wtk.*hess(:,kr2,d1))*h0011_1;
            shh0020(:,d1) = (wtk.*hess(:,kr2,d1))*h0020_1;
            
            for d2=lds.phases
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
                stensv11(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v11_1;
                stensv2(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v2_1;                
                stensv22(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v22_1;                
                for d3=lds.phases
                    stens4v1(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v1_1;
                    stens4v11(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v11_1;
                    stens4v2(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v2_1;
                end
                stens4v1v1(:,d2,d1) = stens4v1v1(:,d2,d1)+(wtk.*stens4v1(:,kr2,d2,d1))*v1_1;
                stens4v1v11(:,d2,d1) = stens4v1v11(:,d2,d1)+(wtk.*stens4v1(:,kr2,d2,d1))*v11_1;
                stens4v1v2(:,d2,d1) = stens4v1v2(:,d2,d1)+(wtk.*stens4v1(:,kr2,d2,d1))*v2_1;
                stens4v2v2(:,d2,d1) = stens4v2v2(:,d2,d1)+(wtk.*stens4v2(:,kr2,d2,d1))*v2_1;
                stens4v11v2(:,d2,d1) = stens4v11v2(:,d2,d1)+(wtk.*stens4v11(:,kr2,d2,d1))*v2_1;
            end    
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
            stenv1v11(:,d1)=stenv1v11(:,d1)+wtk.*stensv1(:,kr2,d1)*v11_1;
            stenv1v2(:,d1)=stenv1v2(:,d1)+wtk.*stensv1(:,kr2,d1)*v2_1;
            stenv1v22(:,d1)=stenv1v22(:,d1)+wtk.*stensv1(:,kr2,d1)*v22_1;
            stenv11v2(:,d1)=stenv11v2(:,d1)+wtk.*stensv11(:,kr2,d1)*v2_1;
            stenv11v22(:,d1)=stenv11v22(:,d1)+wtk.*stensv11(:,kr2,d1)*v22_1;
            stenv2v2(:,d1)=stenv2v2(:,d1)+wtk.*stensv2(:,kr2,d1)*v2_1;              
            stenv2v22(:,d1)=stenv2v22(:,d1)+wtk.*stensv2(:,kr2,d1)*v22_1;              
            stenv22v22(:,d1)=stenv22v22(:,d1)+wtk.*stensv22(:,kr2,d1)*v22_1;              
            
            stens4v1v1v1(:,d1) = stens4v1v1v1(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v1_1;
            stens4v1v1v11(:,d1) = stens4v1v1v11(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v11_1;
            stens4v1v1v2(:,d1) = stens4v1v1v2(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v2_1;
            stens4v1v1v22(:,d1) = stens4v1v1v22(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v22_1;
            stens4v1v11v2(:,d1) = stens4v1v11v2(:,d1)+(wtk.*stens4v1v11(:,kr2,d1))*v2_1;
            stens4v2v2v2(:,d1) = stens4v2v2v2(:,d1)+(wtk.*stens4v2v2(:,kr2,d1))*v2_1;
            stens4v2v2v22(:,d1) = stens4v2v2v22(:,d1)+(wtk.*stens4v2v2(:,kr2,d1))*v22_1;
            stens4v1v2v2(:,d1) = stens4v1v2v2(:,d1)+(wtk.*stens4v1v2(:,kr2,d1))*v2_1;
            stens4v11v2v2(:,d1) = stens4v11v2v2(:,d1)+(wtk.*stens4v11v2(:,kr2,d1))*v2_1;

        end
        fxhessv1h3000(range2)=wtk.*shv1(:,kr2)*h3000_1;        
        fxhessv1h2100(range2)=wtk.*shv1(:,kr2)*h2100_1;
        fxhessv1h1110(range2)=wtk.*shv1(:,kr2)*h1110_1;
        fxhessv1h2010(range2)=wtk.*shv1(:,kr2)*h2010_1;
        fxhessv1h2001(range2)=wtk.*shv1(:,kr2)*h2001_1;        
        fxhessv1h1200(range2)=wtk.*shv1(:,kr2)*h1200_1;        
        fxhessv1h1020(range2)=wtk.*shv1(:,kr2)*h1020_1;        
        fxhessv1h1002(range2)=wtk.*shv1(:,kr2)*h1002_1;
        fxhessv1h1011(range2)=wtk.*shv1(:,kr2)*h1011_1;
        fxhessv1h0111(range2)=wtk.*shv1(:,kr2)*h0111_1;
        fxhessv1h1011(range2)=wtk.*shv1(:,kr2)*h1011_1;
        fxhessv1h1101(range2)=wtk.*shv1(:,kr2)*h1101_1;
        fxhessv1h0021(range2)=wtk.*shv1(:,kr2)*h0021_1;
        fxhessv1h1020(range2)=wtk.*shv1(:,kr2)*h1020_1;
        fxhessv1h0120(range2)=wtk.*shv1(:,kr2)*h0120_1;
        fxhessv11h1020(range2)=wtk.*shv11(:,kr2)*h1020_1;        
        fxhessv11h0021(range2)=wtk.*shv11(:,kr2)*h0021_1;        
        fxhessv11h2001(range2)=wtk.*shv11(:,kr2)*h2001_1;
        fxhessv11h2010(range2)=wtk.*shv11(:,kr2)*h2010_1;        
        fxhessv11h3000(range2)=wtk.*shv11(:,kr2)*h3000_1;
        fxhessv2h3000(range2)=wtk.*shv2(:,kr2)*h3000_1;
        fxhessv2h2100(range2)=wtk.*shv2(:,kr2)*h2100_1;
        fxhessv2h1101(range2)=wtk.*shv2(:,kr2)*h1101_1;            
        fxhessv2h0012(range2)=wtk.*shv2(:,kr2)*h0012_1;
        fxhessv2h2001(range2)=wtk.*shv2(:,kr2)*h2001_1;
        fxhessv2h2010(range2)=wtk.*shv2(:,kr2)*h2010_1;        
        fxhessv2h1011(range2)=wtk.*shv2(:,kr2)*h1011_1;
        fxhessv2h0111(range2)=wtk.*shv2(:,kr2)*h0111_1;
        fxhessv2h0021(range2)=wtk.*shv2(:,kr2)*h0021_1;        
        fxhessv2h1110(range2)=wtk.*shv2(:,kr2)*h1110_1;        
        fxhessv22h0120(range2)=wtk.*shv22(:,kr2)*h0120_1;
        fxhessv22h0030(range2)=wtk.*shv22(:,kr2)*h0030_1;
        fxhessv22h3000(range2)=wtk.*shv22(:,kr2)*h3000_1;
        fxhessv22h2001(range2)=wtk.*shv22(:,kr2)*h2001_1;
        fxhessv22h2100(range2)=wtk.*shv22(:,kr2)*h2100_1;
        fxhessv22h2010(range2)=wtk.*shv22(:,kr2)*h2010_1;
        fxhessv22h1020(range2)=wtk.*shv22(:,kr2)*h1020_1;

        fxhessh0011h0011(range2)=wtk.*shh0011(:,kr2)*h0011_1;
        fxhessh0011h1100(range2)=wtk.*shh0011(:,kr2)*h1100_1;
        fxhessh0020h0011(range2)=wtk.*shh0020(:,kr2)*h0011_1;
        fxhessh0020h1001(range2)=wtk.*shh0020(:,kr2)*h1001_1;       
        fxhessh0020h0002(range2)=wtk.*shh0020(:,kr2)*h0002_1;
        fxhessh0020h0101(range2)=wtk.*shh0020(:,kr2)*h0101_1;
        fxhessh0020h1100(range2)=wtk.*shh0020(:,kr2)*h1100_1;        
        fxhessh0101h1010(range2)=wtk.*shh0101(:,kr2)*h1010_1;
        fxhessh0110h0011(range2)=wtk.*shh0110(:,kr2)*h0011_1;
        fxhessh0110h1010(range2)=wtk.*shh0110(:,kr2)*h1010_1;
        fxhessh0110h1001(range2)=wtk.*shh0110(:,kr2)*h1001_1;
        fxhessh1001h1001(range2)=wtk.*shh1001(:,kr2)*h1001_1;
        fxhessh1001h1100(range2)=wtk.*shh1001(:,kr2)*h1100_1;
        fxhessh1001h1010(range2)=wtk.*shh1001(:,kr2)*h1010_1;
        fxhessh1010h0011(range2)=wtk.*shh1010(:,kr2)*h0011_1;       
        fxhessh1010h1100(range2)=wtk.*shh1010(:,kr2)*h1100_1;
        fxhessh1010h1010(range2)=wtk.*shh1010(:,kr2)*h1010_1;
        fxhessh1100h1100(range2)=wtk.*shh1100(:,kr2)*h1100_1;
        
        fxhessh2000h0110(range2)=wtk.*shh2000(:,kr2)*h0110_1;        
        fxhessh2000h1100(range2)=wtk.*shh2000(:,kr2)*h1100_1;        
        fxhessh2000h1010(range2)=wtk.*shh2000(:,kr2)*h1010_1;        
        fxhessh2000h1001(range2)=wtk.*shh2000(:,kr2)*h1001_1;        
        fxhessh2000h0200(range2)=wtk.*shh2000(:,kr2)*h0200_1;
        fxhessh2000h0002(range2)=wtk.*shh2000(:,kr2)*h0002_1;
        fxhessh2000h0101(range2)=wtk.*shh2000(:,kr2)*h0101_1;
        fxhessh2000h0011(range2)=wtk.*shh2000(:,kr2)*h0011_1;
        fxhessh2000h2000(range2)=wtk.*shh2000(:,kr2)*h2000_1;
        fxhessh2000h0020(range2)=wtk.*shh2000(:,kr2)*h0020_1;
              
        fxtensv1v1h1001(range2)=wtk.*stenv1v1(:,kr2)*h1001_1;
        fxtensv1v1h0200(range2)=wtk.*stenv1v1(:,kr2)*h0200_1;
        fxtensv1v1h0020(range2)=wtk.*stenv1v1(:,kr2)*h0020_1;
        fxtensv1v1h0002(range2)=wtk.*stenv1v1(:,kr2)*h0002_1;                
        fxtensv1v1h0110(range2)=wtk.*stenv1v1(:,kr2)*h0110_1;
        fxtensv1v1h0101(range2)=wtk.*stenv1v1(:,kr2)*h0101_1;        
        fxtensv1v1h0011(range2)=wtk.*stenv1v1(:,kr2)*h0011_1;      
        fxtensv1v1h2000(range2)=wtk.*stenv1v1(:,kr2)*h2000_1;
        fxtensv1v1h1100(range2)=wtk.*stenv1v1(:,kr2)*h1100_1;
        fxtensv1v1h1010(range2)=wtk.*stenv1v1(:,kr2)*h1010_1;
        fxtensv1v11h2000(range2)=wtk.*stenv1v11(:,kr2)*h2000_1;        
        fxtensv1v11h0011(range2)=wtk.*stenv1v11(:,kr2)*h0011_1;        
        fxtensv1v11h0020(range2)=wtk.*stenv1v11(:,kr2)*h0020_1;             
        fxtensv1v11h1001(range2)=wtk.*stenv1v11(:,kr2)*h1001_1;        
        fxtensv1v11h1010(range2)=wtk.*stenv1v11(:,kr2)*h1010_1;
        fxtensv1v11h1100(range2)=wtk.*stenv1v11(:,kr2)*h1100_1;

        fxtensv1v2h2000(range2)=wtk.*stenv1v2(:,kr2)*h2000_1;
        fxtensv1v2h1100(range2)=wtk.*stenv1v2(:,kr2)*h1100_1;        
        fxtensv1v2h1001(range2)=wtk.*stenv1v2(:,kr2)*h1001_1;        
        fxtensv1v2h1010(range2)=wtk.*stenv1v2(:,kr2)*h1010_1;                        
        fxtensv1v2h0101(range2)=wtk.*stenv1v2(:,kr2)*h0101_1;        
        fxtensv1v2h0110(range2)=wtk.*stenv1v2(:,kr2)*h0110_1;        
        fxtensv1v2h0011(range2)=wtk.*stenv1v2(:,kr2)*h0011_1;        
        fxtensv1v22h0020(range2)=wtk.*stenv1v22(:,kr2)*h0020_1;                
        fxtensv1v22h0110(range2)=wtk.*stenv1v22(:,kr2)*h0110_1;                
        fxtensv1v22h1001(range2)=wtk.*stenv1v22(:,kr2)*h1001_1;                
        fxtensv1v22h1100(range2)=wtk.*stenv1v22(:,kr2)*h1100_1;        
        fxtensv1v22h1010(range2)=wtk.*stenv1v22(:,kr2)*h1010_1;        
        fxtensv1v22h2000(range2)=wtk.*stenv1v22(:,kr2)*h2000_1;
        fxtensv11v2h2000(range2)=wtk.*stenv11v2(:,kr2)*h2000_1;
        fxtensv11v2h0011(range2)=wtk.*stenv11v2(:,kr2)*h0011_1;        
        fxtensv11v2h1010(range2)=wtk.*stenv11v2(:,kr2)*h1010_1;                
        fxtensv11v22h0020(range2)=wtk.*stenv11v22(:,kr2)*h0020_1;        
        fxtensv11v22h2000(range2)=wtk.*stenv11v22(:,kr2)*h2000_1;        
        fxtensv2v2h1100(range2)=wtk.*stenv2v2(:,kr2)*h1100_1;       
        fxtensv2v2h2000(range2)=wtk.*stenv2v2(:,kr2)*h2000_1;        
        fxtensv2v2h0101(range2)=wtk.*stenv2v2(:,kr2)*h0101_1;           
        fxtensv2v2h1001(range2)=wtk.*stenv2v2(:,kr2)*h1001_1;                
        fxtensv2v2h0002(range2)=wtk.*stenv2v2(:,kr2)*h0002_1;        
        fxtensv2v2h0011(range2)=wtk.*stenv2v2(:,kr2)*h0011_1;        
        fxtensv2v22h1100(range2)=wtk.*stenv2v22(:,kr2)*h1100_1;        
        fxtensv2v22h0020(range2)=wtk.*stenv2v22(:,kr2)*h0020_1;        
        fxtensv2v22h2000(range2)=wtk.*stenv2v22(:,kr2)*h2000_1;                
        fxtensv2v22h0011(range2)=wtk.*stenv2v22(:,kr2)*h0011_1;        
        fxtensv2v22h1010(range2)=wtk.*stenv2v22(:,kr2)*h1010_1;        
        fxtensv2v22h0110(range2)=wtk.*stenv2v22(:,kr2)*h0110_1;        
        fxtensv22v22h2000(range2)=wtk.*stenv22v22(:,kr2)*h2000_1;                
                   
        fxtens4v1v1v1v1(range2)=wtk.*stens4v1v1v1(:,kr2)*v1_1;
        fxtens4v1v1v1v11(range2)=wtk.*stens4v1v1v1(:,kr2)*v11_1;
        fxtens4v1v1v1v2(range2)=wtk.*stens4v1v1v1(:,kr2)*v2_1;
        fxtens4v1v1v1v22(range2)=wtk.*stens4v1v1v1(:,kr2)*v22_1;
        fxtens4v1v1v11v11(range2)=wtk.*stens4v1v1v11(:,kr2)*v11_1;        
        fxtens4v1v1v11v22(range2)=wtk.*stens4v1v1v11(:,kr2)*v22_1;        
        fxtens4v1v1v11v2(range2)=wtk.*stens4v1v1v11(:,kr2)*v2_1;        
        fxtens4v1v1v2v2(range2)=wtk.*stens4v1v1v2(:,kr2)*v2_1;        
        fxtens4v1v1v2v22(range2)=wtk.*stens4v1v1v2(:,kr2)*v22_1;        
        fxtens4v1v1v22v22(range2)=wtk.*stens4v1v1v22(:,kr2)*v22_1;        
        fxtens4v1v1v22v22(range2)=wtk.*stens4v1v1v22(:,kr2)*v22_1;        
        fxtens4v1v11v2v2(range2)=wtk.*stens4v1v11v2(:,kr2)*v2_1;       
        fxtens4v1v11v2v22(range2)=wtk.*stens4v1v11v2(:,kr2)*v22_1;        
        fxtens4v1v2v2v22(range2)=wtk.*stens4v1v2v2(:,kr2)*v22_1;               
        fxtens4v11v2v2v22(range2)=wtk.*stens4v11v2v2(:,kr2)*v22_1;                
        fxtens4v2v2v2v22(range2)=wtk.*stens4v2v2v2(:,kr2)*v22_1;        
        fxtens4v2v2v22v22(range2)=wtk.*stens4v2v2v22(:,kr2)*v22_1;                
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end

h2000C = reshape(h2000Cps,lds.ncoords-lds.nphase,1);
h0020C = reshape(h0020Cps,lds.ncoords-lds.nphase,1);
h1100C = reshape(h1100Cps,lds.ncoords-lds.nphase,1);
h1010C = reshape(h1010Cps,lds.ncoords-lds.nphase,1);
h1011C = reshape(h1011Cps,lds.ncoords-lds.nphase,1);
h1001C = reshape(h1001Cps,lds.ncoords-lds.nphase,1);
h0011C = reshape(h0011Cps,lds.ncoords-lds.nphase,1);
h0110C = conj(h1001C);
h2100C = reshape(h2100Cps,lds.ncoords-lds.nphase,1);
h0021C = reshape(h0021Cps,lds.ncoords-lds.nphase,1);
h1110C = reshape(h1110Cps,lds.ncoords-lds.nphase,1);

%computation of h4000
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v1v1.'+6*T*fxtensv1v1h2000.'+3*T*fxhessh2000h2000.'+4*T*fxhessv1h3000.'; rl(lds.ncoords) = 0;
h4000 = J4theta1plus\rl.';

%computation of h3100
diffh2000 = T*fxjac*h2000-2*i*theta1*h2000C+diffh2000.';
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v1v11.'+3*T*fxtensv1v1h1100.'+3*T*fxtensv1v11h2000.'+T*fxhessv11h3000.'...
    +3*T*fxhessv1h2100.'+3*T*fxhessh2000h1100.'-3*T*alpha11001*diffh2000-6*T*a21001*h2000C; rl(lds.ncoords) = 0;
h3100 = J2theta1plus\rl.';

%computation of h3100
diffh0020 = T*fxjac*h0020-2*i*theta2*h0020C+diffh0020.';
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v2v2v2v22.'+3*T*fxtensv2v2h0011.'+3*T*fxtensv2v22h0020.'+T*fxhessv22h0030.'...
    +3*T*fxhessv2h0021.'+3*T*fxhessh0020h0011.'-3*T*alpha00111*diffh0020-6*T*b00211*h0020C; rl(lds.ncoords) = 0;
h0031 = J2theta2plus\rl.';

%computation of h3010
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v1v2.'+3*T*fxtensv1v1h1010.'+3*T*fxtensv1v2h2000.'+T*fxhessv2h3000.'...
    +3*T*fxhessv1h2010.'+3*T*fxhessh2000h1010.'; rl(lds.ncoords) = 0;
h3010 = J3theta1plustheta2plus\rl.';

%computation of h3001
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v1v22.'+3*T*fxtensv1v1h1001.'+3*T*fxtensv1v22h2000.'+T*fxhessv22h3000.'...
    +3*T*fxhessv1h2001.'+3*T*fxhessh2000h1001.'; rl(lds.ncoords) = 0;
h3001 = J3theta1plustheta2min\rl.';

diffh1100 = T*fxjac*h1100+diffh1100.';
alpha22001= 1/4*phi'*(real(fxtens4v1v1v11v11)'+2*real(fxtensv1v1h0200)'+4*real(fxtensv1v11h1100)'...
    +real(fxhessh2000h0200)'+2*fxhessh1100h1100'+4*real(fxhessv1h1200)'-4*alpha11001*diffh1100);
alpha2200 = alpha22001/T^2;

diffh0011 = T*fxjac*h0011+diffh0011.';
alpha00221= 1/4*phi'*(real(fxtens4v2v2v22v22)'+2*real(fxtensv2v2h0002)'+4*real(fxtensv2v22h0011)'...
    +real(fxhessh0020h0002)'+2*fxhessh0011h0011'+4*real(fxhessv2h0012)'-4*alpha00111*diffh0011);
alpha0022 = alpha00221/T^2;

%computation of h2200
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxtens4v1v1v11v11)'+2*T*real(fxtensv1v1h0200)'+4*T*real(fxtensv1v11h1100)'...
    +T*real(fxhessh2000h0200)'+2*T*fxhessh1100h1100'+4*T*real(fxhessv1h1200)'-8*T*real(a21001)*h1100C-4*T*alpha22001*f1...
    -4*T*alpha11001*diffh1100; rl(lds.ncoords+1) = 0;
h2200 = Jjac3\rl.';
h2200 = h2200(lds.coords);

%computation of h0022
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxtens4v2v2v22v22)'+2*T*real(fxtensv2v2h0002)'+4*T*real(fxtensv2v22h0011)'...
    +T*real(fxhessh0020h0002)'+2*T*fxhessh0011h0011'+4*T*real(fxhessv2h0012)'-8*T*real(b00211)*h0011C-4*T*alpha00221*f1...
    -4*T*alpha00111*diffh0011; rl(lds.ncoords+1) = 0;
h0022 = Jjac3\rl.';
h0022 = h0022(lds.coords);

%computation of h2020
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v2v2.'+T*fxtensv1v1h0020.'+T*fxtensv2v2h2000.'+4*T*fxtensv1v2h1010.'...
    +T*fxhessh2000h0020.'+2*T*fxhessv2h2010.'+2*T*fxhessh1010h1010.'+2*T*fxhessv1h1020.'; rl(lds.ncoords) = 0;
h2020 = J2theta1plus2theta2plus\rl.';

%computation of h2002
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v22v22.'+T*fxtensv22v22h2000.'+4*T*fxtensv1v22h1001.'+T*fxtensv1v1h0002.'...
    +2*T*fxhessv22h2001.'+T*fxhessh2000h0002.'+2*T*fxhessh1001h1001.'+2*T*fxhessv1h1002.'; rl(lds.ncoords) = 0;
h2002 = J2theta1plus2theta2min\rl.';

%computation of h2110
diffh1010 = T*fxjac*h1010-i*theta1*h1010C-i*theta2*h1010C+diffh1010.';
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v11v2.'+T*fxtensv1v1h0110.'+2*T*fxtensv1v11h1010.'+T*fxtensv11v2h2000.'...
    +2*T*fxtensv1v2h1100.'+T*fxhessv11h2010.'+2*T*fxhessh1010h1100.'+T*fxhessv2h2100.'+T*fxhessh2000h0110.'...
    +2*T*fxhessv1h1110.'-2*T*a21001*h1010C-2*T*b11101*h1010C-2*T*alpha11001*diffh1010; rl(lds.ncoords) = 0;
h2110 = Jtheta1plustheta2plus\rl.';

%computation of h1021
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v2v2v22.'+T*fxtensv2v2h1001.'+2*T*fxtensv2v22h1010.'+T*fxtensv1v22h0020.'...
    +2*T*fxtensv1v2h0011.'+T*fxhessv22h1020.'+2*T*fxhessh1010h0011.'+T*fxhessv1h0021.'+T*fxhessh0020h1001.'...
    +2*T*fxhessv2h1011.'-2*T*b00211*h1010C-2*T*a10111*h1010C-2*T*alpha00111*diffh1010; rl(lds.ncoords) = 0;
h1021 = Jtheta1plustheta2plus\rl.';

%computation of h2101
diffh1001 = T*fxjac*h1001-i*theta1*h1001C+i*theta2*h1001C+diffh1001.';
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v11v22.'+T*fxtensv1v1h0101.'+2*T*fxtensv1v11h1001.'+T*fxtensv11v22h2000.'...
    +2*T*fxtensv1v22h1100.'+2*T*fxhessh1001h1100.'+2*T*fxhessv1h1101.'+T*fxhessv22h2100.'+T*fxhessv11h2001.'...
    +T*fxhessh2000h0101.'-2*T*a21001*h1001C-2*T*b11011*h1001C-2*T*alpha11001*diffh1001; rl(lds.ncoords) = 0;
h2101 = Jtheta1plustheta2min\rl.';

%computation of h0121
diffh0110 = conj(diffh1001);
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v11v2v2v22.'+T*fxtensv2v2h0101.'+2*T*fxtensv2v22h0110.'+T*fxtensv11v22h0020.'...
    +2*T*fxtensv11v2h0011.'+2*T*fxhessh0110h0011.'+2*T*fxhessv2h0111.'+T*fxhessv11h0021.'+T*fxhessv22h0120.'...
    +T*fxhessh0020h0101.'-2*T*b00211*h0110C-2*T*b01111*h0110C-2*T*alpha00111*diffh0110; rl(lds.ncoords) = 0;
h0121 = Jtheta1mintheta2plus\rl.';

%computation of h2011
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v2v22.'+T*fxtensv1v1h0011.'+2*T*fxtensv1v2h1001.'+T*fxtensv2v22h2000.'...
    +2*T*fxtensv1v22h1010.'+T*fxhessv22h2010.'+T*fxhessv2h2001.'+T*fxhessh2000h0011.'+2*T*fxhessh1001h1010.'...
    +2*T*fxhessv1h1011.'-2*T*a10111*h2000C-T*alpha00111*diffh2000; rl(lds.ncoords) = 0;
h2011 = J2theta1plus\rl.';

%computation of h1120
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v11v2v2.'+T*fxtensv2v2h1100.'+2*T*fxtensv1v2h0110.'+T*fxtensv1v11h0020.'...
    +2*T*fxtensv11v2h1010.'+T*fxhessv11h1020.'+T*fxhessv1h0120.'+T*fxhessh0020h1100.'+2*T*fxhessh0110h1010.'...
    +2*T*fxhessv2h1110.'-2*T*b11101*h0020C-T*alpha11001*diffh0020; rl(lds.ncoords) = 0;
h1120 = J2theta2plus\rl.';

%computation of alpha11111
alpha11111 = phi'*(real(fxtens4v1v11v2v22)'+real(fxtensv1v11h0011)'+2*real(fxtensv1v2h0101)'+2*real(fxtensv1v22h0110)'...
    +real(fxtensv2v22h1100)'+2*real(fxhessv1h0111)'+real(fxhessh0110h1001)'+real(fxhessh0101h1010)'+fxhessh0011h1100.'...
    +2*real(fxhessv2h1101)'-alpha00111*diffh1100-alpha11001*diffh0011);
alpha1111 = alpha11111/T^2;

%computation of h1111
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxtens4v1v11v2v22)'+T*real(fxtensv1v11h0011)'+2*T*real(fxtensv1v2h0101)'...
    +2*T*real(fxtensv1v22h0110)'+T*real(fxtensv2v22h1100)'+2*T*real(fxhessv1h0111)'+T*real(fxhessh0110h1001)'...
    +T*real(fxhessh0101h1010)'+T*fxhessh0011h1100'+2*T*real(fxhessv2h1101)'-2*T*real(a01111)*h1100C-T*alpha00111*diffh1100...
    -2*T*real(b11011)*h0011C-T*alpha11111*f1-T*alpha11001*diffh0011; rl(lds.ncoords+1) = 0;
h1111 = Jjac3\rl.';
h1111 = h1111(lds.coords);


% function
range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    
    v1_1 = v1(range4);
    v11_1 = conj(v1_1);
    v2_1 = v2(range4);
    v22_1 = conj(v2_1);
    h0020_1 = h0020(range4);
    h0002_1 = conj(h0020_1);
    h0011_1 = h0011(range4);
    h0021_1 = h0021(range4);
    h0012_1 = conj(h0021_1);        
    h0022_1 = h0022(range4);
    h0030_1 = h0030(range4);
    h0031_1 = h0031(range4);
    h1010_1 = h1010(range4);        
    h1001_1 = h1001(range4);
    h1011_1 = h1011(range4);        
    h0101_1 = conj(h1010_1);    
    h0110_1 = conj(h1001_1);
    h0111_1 = h0111(range4);    
    h0121_1 = h0121(range4);
    h2000_1 = h2000(range4);    
    h0200_1 = conj(h2000_1);            
    h1020_1 = h1020(range4);      
    h0120_1 = h0120(range4);        
    h1002_1 = conj(h0120_1);        
    h1021_1 = h1021(range4);
    h1012_1 = conj(h0121_1);        
    h1100_1 = h1100(range4);        
    h1110_1 = h1110(range4);    
    h1101_1 = h1101(range4);            
    h1111_1 = h1111(range4);
    h1120_1 = h1120(range4);
    h2100_1 = h2100(range4);                
    h1200_1 = conj(h2100_1);    
    h2101_1 = h2101(range4);
    h1210_1 = conj(h2101_1);
    h2010_1 = h2010(range4);     
    h2001_1 = h2001(range4);    
    h0210_1 = conj(h2001_1);
    h2011_1 = h2011(range4);        
    h2110_1 = h2110(range4);
    h2200_1 = h2200(range4);
    h3000_1 = h3000(range4);
    h3100_1 = h3100(range4);
    
    range5 = lds.phases;         
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        stens5v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase,lds.nphase);
        stens5v1v11 = zeros(lds.nphase,lds.nphase,lds.nphase,lds.nphase);
        stens5v1v2 = zeros(lds.nphase,lds.nphase,lds.nphase,lds.nphase);
        stens5v2v2 = zeros(lds.nphase,lds.nphase,lds.nphase,lds.nphase);
        
        stens4v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v1v11 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v2v2 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v2v22 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v1v22 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v1v2 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v11v11 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v11v2 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens5v1v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens5v1v2v2 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens5v1v1v11 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens5v2v2v2 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens5v1v11v2 = zeros(lds.nphase,lds.nphase,lds.nphase);
        
        stenv1v1 = zeros(lds.nphase,lds.nphase);
        stenv1v11 = zeros(lds.nphase,lds.nphase);
        stenv1v2 = zeros(lds.nphase,lds.nphase);
        stenv1v22 = zeros(lds.nphase,lds.nphase);
        stenv11v11 = zeros(lds.nphase,lds.nphase);
        stenv11h2000 = zeros(lds.nphase,lds.nphase);
        stenv1h0200 = zeros(lds.nphase,lds.nphase);
        stenv2v2 = zeros(lds.nphase,lds.nphase);
        stenv22v22 = zeros(lds.nphase,lds.nphase);
        stenv22h0020 = zeros(lds.nphase,lds.nphase);
        stenv2h0002 = zeros(lds.nphase,lds.nphase);
        stenv1h0020 = zeros(lds.nphase,lds.nphase);
        stenv22h1001 = zeros(lds.nphase,lds.nphase);
        stenv1h0011 = zeros(lds.nphase,lds.nphase);
        stenv1h0110 = zeros(lds.nphase,lds.nphase);
        stenv1h1010 = zeros(lds.nphase,lds.nphase);
        stenv11h1010 = zeros(lds.nphase,lds.nphase);
        stenv2h1100 = zeros(lds.nphase,lds.nphase);
        stenv1h1001 = zeros(lds.nphase,lds.nphase);
        stenv2v22 = zeros(lds.nphase,lds.nphase);
        stenv1h0022 = zeros(lds.nphase,lds.nphase);
        stenv1h1100 = zeros(lds.nphase,lds.nphase);
        stenv11v22 = zeros(lds.nphase,lds.nphase);
        stenv11v2 = zeros(lds.nphase,lds.nphase);
        stenv22h2000 = zeros(lds.nphase,lds.nphase);
        stenv22h1010 = zeros(lds.nphase,lds.nphase);
        stenv2h2000 = zeros(lds.nphase,lds.nphase);
        stenv11h1001 = zeros(lds.nphase,lds.nphase);
        stenv2h0110 = zeros(lds.nphase,lds.nphase);
        stenv2h1001 = zeros(lds.nphase,lds.nphase);
        stenv2h0011 = zeros(lds.nphase,lds.nphase);
        stenv2h1010 = zeros(lds.nphase,lds.nphase);
        stenv22h0110 = zeros(lds.nphase,lds.nphase);
        stenv11h0020 = zeros(lds.nphase,lds.nphase);
        stenv11h0110 = zeros(lds.nphase,lds.nphase);
        stens4v1v1v1 = zeros(lds.nphase,lds.nphase);
        stens4v1v1v2 = zeros(lds.nphase,lds.nphase);
        stens4v1v1v22 = zeros(lds.nphase,lds.nphase);
        stens4v1v11v2 = zeros(lds.nphase,lds.nphase);
        stens4v1v11v11 = zeros(lds.nphase,lds.nphase);
        stens4v1v1v11 = zeros(lds.nphase,lds.nphase);
        stens4v2v2v2 = zeros(lds.nphase,lds.nphase);
        stens4v2v22v22 = zeros(lds.nphase,lds.nphase);
        stens4v2v2v22 = zeros(lds.nphase,lds.nphase);
        stens4v1v22v22 = zeros(lds.nphase,lds.nphase);
        stens4v1v2v2 = zeros(lds.nphase,lds.nphase);
        stens4v1v2v22 = zeros(lds.nphase,lds.nphase);
        stens4v11v11v2 = zeros(lds.nphase,lds.nphase);
        stens4v1v11v22 = zeros(lds.nphase,lds.nphase);
        stens4v11v2v22 = zeros(lds.nphase,lds.nphase);
        stens4v11v2v2 = zeros(lds.nphase,lds.nphase);
        stens5v1v1v11v11 = zeros(lds.nphase,lds.nphase);
        stens5v2v2v2v22 = zeros(lds.nphase,lds.nphase);
        stens5v1v2v2v22 = zeros(lds.nphase,lds.nphase);
        stens5v1v1v11v2 = zeros(lds.nphase,lds.nphase);
        stens5v1v11v2v2 = zeros(lds.nphase,lds.nphase);
        stens5v1v1v1v11 = zeros(lds.nphase,lds.nphase);
        
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        tens4 = ctens4(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,xt,pt,lds.ActiveParams);
        tens5 = ctens5(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,lds.Der5,xt,pt,lds.ActiveParams);
        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shv11(:,d1) = (wtk.*hess(:,kr2,d1))*v11_1;
            shv2(:,d1) = (wtk.*hess(:,kr2,d1))*v2_1;
            shv22(:,d1) = (wtk.*hess(:,kr2,d1))*v22_1;
            shh2000(:,d1) = (wtk.*hess(:,kr2,d1))*h2000_1;
            shh0200(:,d1) = (wtk.*hess(:,kr2,d1))*h0200_1;
            shh2100(:,d1) = (wtk.*hess(:,kr2,d1))*h2100_1;
            shh0021(:,d1) = (wtk.*hess(:,kr2,d1))*h0021_1;
            shh0020(:,d1) = (wtk.*hess(:,kr2,d1))*h0020_1;
            shh0012(:,d1) = (wtk.*hess(:,kr2,d1))*h0012_1;
            shh0011(:,d1) = (wtk.*hess(:,kr2,d1))*h0011_1;
            shh1200(:,d1) = (wtk.*hess(:,kr2,d1))*h1200_1;
            shh1100(:,d1) = (wtk.*hess(:,kr2,d1))*h1100_1;
            shh0101(:,d1) = (wtk.*hess(:,kr2,d1))*h0101_1;
            shh2001(:,d1) = (wtk.*hess(:,kr2,d1))*h2001_1;
            shh1011(:,d1) = (wtk.*hess(:,kr2,d1))*h1011_1;
            shh1010(:,d1) = (wtk.*hess(:,kr2,d1))*h1010_1;
            shh1001(:,d1) = (wtk.*hess(:,kr2,d1))*h1001_1;
            shh0120(:,d1) = (wtk.*hess(:,kr2,d1))*h0120_1;
            shh1110(:,d1) = (wtk.*hess(:,kr2,d1))*h1110_1;
            shh0110(:,d1) = (wtk.*hess(:,kr2,d1))*h0110_1;
            shh0002(:,d1) = (wtk.*hess(:,kr2,d1))*h0002_1;
            
            for d2=lds.phases
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
                stensv11(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v11_1;
                stensv2(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v2_1;                
                stensv22(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v22_1;                
                for d3=lds.phases
                    stens4v1(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v1_1;
                    stens4v2(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v2_1;
                    stens4v11(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v11_1;
                    for d4=lds.phases
                        stens5v1(:,d4,d3,d2,d1) = (wtk.*tens5(:,kr2,d1,d2,d3,d4))*v1_1;
                        stens5v2(:,d4,d3,d2,d1) = (wtk.*tens5(:,kr2,d1,d2,d3,d4))*v2_1;
                    end
                    stens5v1v1(:,d3,d2,d1) = stens5v1v1(:,d3,d2,d1)+wtk.*stens5v1(:,kr2,d3,d2,d1)*v1_1;
                    stens5v1v11(:,d3,d2,d1) = stens5v1v11(:,d3,d2,d1)+wtk.*stens5v1(:,kr2,d3,d2,d1)*v11_1;
                    stens5v1v2(:,d3,d2,d1) = stens5v1v2(:,d3,d2,d1)+wtk.*stens5v1(:,kr2,d3,d2,d1)*v2_1;
                    stens5v2v2(:,d3,d2,d1) = stens5v2v2(:,d3,d2,d1)+wtk.*stens5v2(:,kr2,d3,d2,d1)*v2_1;
                end
                stens4v1v1(:,d2,d1) = stens4v1v1(:,d2,d1)+(wtk.*stens4v1(:,kr2,d2,d1))*v1_1;
                stens4v1v11(:,d2,d1) = stens4v1v11(:,d2,d1)+(wtk.*stens4v1(:,kr2,d2,d1))*v11_1;
                stens4v2v2(:,d2,d1) = stens4v2v2(:,d2,d1)+(wtk.*stens4v2(:,kr2,d2,d1))*v2_1;
                stens4v2v22(:,d2,d1) = stens4v2v22(:,d2,d1)+(wtk.*stens4v2(:,kr2,d2,d1))*v22_1;
                stens4v1v22(:,d2,d1) = stens4v1v22(:,d2,d1)+(wtk.*stens4v1(:,kr2,d2,d1))*v22_1;
                stens4v1v2(:,d2,d1) = stens4v1v2(:,d2,d1)+(wtk.*stens4v1(:,kr2,d2,d1))*v2_1;
                stens4v11v11(:,d2,d1) = stens4v11v11(:,d2,d1)+(wtk.*stens4v11(:,kr2,d2,d1))*v11_1;
                stens4v11v2(:,d2,d1) = stens4v11v2(:,d2,d1)+(wtk.*stens4v11(:,kr2,d2,d1))*v2_1;
                
                stens5v1v1v1(:,d2,d1) = stens5v1v1v1(:,d2,d1)+wtk.*stens5v1v1(:,kr2,d2,d1)*v1_1;
                stens5v1v2v2(:,d2,d1) = stens5v1v2v2(:,d2,d1)+wtk.*stens5v1v2(:,kr2,d2,d1)*v2_1;
                stens5v1v1v11(:,d2,d1) = stens5v1v1v11(:,d2,d1)+wtk.*stens5v1v1(:,kr2,d2,d1)*v11_1;
                stens5v2v2v2(:,d2,d1) = stens5v2v2v2(:,d2,d1)+wtk.*stens5v2v2(:,kr2,d2,d1)*v2_1;
                stens5v1v11v2(:,d2,d1) = stens5v1v11v2(:,d2,d1)+wtk.*stens5v1v11(:,kr2,d2,d1)*v2_1;

            end

            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
            stenv1v11(:,d1)=stenv1v11(:,d1)+wtk.*stensv1(:,kr2,d1)*v11_1;
            stenv1v2(:,d1)=stenv1v2(:,d1)+wtk.*stensv1(:,kr2,d1)*v2_1;
            stenv1v22(:,d1)=stenv1v22(:,d1)+wtk.*stensv1(:,kr2,d1)*v22_1;
            stenv11v11(:,d1)=stenv11v11(:,d1)+wtk.*stensv11(:,kr2,d1)*v11_1;
            stenv11h2000(:,d1)=stenv11h2000(:,d1)+wtk.*stensv11(:,kr2,d1)*h2000_1;
            stenv1h0200(:,d1)=stenv1h0200(:,d1)+wtk.*stensv1(:,kr2,d1)*h0200_1;              
            stenv2h0011(:,d1)=stenv2h0011(:,d1)+wtk.*stensv2(:,kr2,d1)*h0011_1;              
            stenv2v2(:,d1)=stenv2v2(:,d1)+wtk.*stensv2(:,kr2,d1)*v2_1; 
            stenv22v22(:,d1)=stenv22v22(:,d1)+wtk.*stensv22(:,kr2,d1)*v22_1;
            stenv22h0020(:,d1)=stenv22h0020(:,d1)+wtk.*stensv22(:,kr2,d1)*h0020_1;
            stenv2h0002(:,d1)=stenv2h0002(:,d1)+wtk.*stensv2(:,kr2,d1)*h0002_1;
            stenv1h0020(:,d1)=stenv1h0020(:,d1)+wtk.*stensv1(:,kr2,d1)*h0020_1;              
            stenv22h1001(:,d1)=stenv22h1001(:,d1)+wtk.*stensv22(:,kr2,d1)*h1001_1;              
            stenv2h1001(:,d1)=stenv2h1001(:,d1)+wtk.*stensv2(:,kr2,d1)*h1001_1; 
            stenv2h1010(:,d1)=stenv2h1010(:,d1)+wtk.*stensv2(:,kr2,d1)*h1010_1; 
            stenv1h0011(:,d1)=stenv1h0011(:,d1)+wtk.*stensv1(:,kr2,d1)*h0011_1;
            stenv1h0110(:,d1)=stenv1h0110(:,d1)+wtk.*stensv1(:,kr2,d1)*h0110_1;
            stenv1h1010(:,d1)=stenv1h1010(:,d1)+wtk.*stensv1(:,kr2,d1)*h1010_1;
            stenv11h1010(:,d1)=stenv11h1010(:,d1)+wtk.*stensv11(:,kr2,d1)*h1010_1;              
            stenv2h1100(:,d1)=stenv2h1100(:,d1)+wtk.*stensv2(:,kr2,d1)*h1100_1;              
            stenv1h1001(:,d1)=stenv1h1001(:,d1)+wtk.*stensv1(:,kr2,d1)*h1001_1;
            stenv2v22(:,d1)=stenv2v22(:,d1)+wtk.*stensv2(:,kr2,d1)*v22_1; 
            stenv2h1001(:,d1)=stenv2h1001(:,d1)+wtk.*stensv2(:,kr2,d1)*h1001_1;
            stenv1h1100(:,d1)=stenv1h1100(:,d1)+wtk.*stensv1(:,kr2,d1)*h1100_1;
            stenv11v22(:,d1)=stenv11v22(:,d1)+wtk.*stensv11(:,kr2,d1)*v22_1;
            stenv11v2(:,d1)=stenv11v2(:,d1)+wtk.*stensv11(:,kr2,d1)*v2_1;              
            stenv22h1010(:,d1)=stenv22h1010(:,d1)+wtk.*stensv22(:,kr2,d1)*h1010_1;              
            stenv2h2000(:,d1)=stenv2h2000(:,d1)+wtk.*stensv2(:,kr2,d1)*h2000_1;
            stenv11h1001(:,d1)=stenv11h1001(:,d1)+wtk.*stensv11(:,kr2,d1)*h1001_1;
            stenv2h0110(:,d1)=stenv2h0110(:,d1)+wtk.*stensv2(:,kr2,d1)*h0110_1; 
            stenv2h1001(:,d1)=stenv2h1001(:,d1)+wtk.*stensv2(:,kr2,d1)*h1001_1;
            stenv2h0011(:,d1)=stenv2h0011(:,d1)+wtk.*stensv2(:,kr2,d1)*h0011_1;
            stenv2h1010(:,d1)=stenv2h1010(:,d1)+wtk.*stensv2(:,kr2,d1)*h1010_1;
            stenv22h0110(:,d1)=stenv22h0110(:,d1)+wtk.*stensv22(:,kr2,d1)*h0110_1;              
            stenv11h0020(:,d1)=stenv11h0020(:,d1)+wtk.*stensv11(:,kr2,d1)*h0020_1;              
            stenv11h0110(:,d1)=stenv11h0110(:,d1)+wtk.*stensv11(:,kr2,d1)*h0110_1;              
            stenv1h0022(:,d1)=stenv1h0022(:,d1)+wtk.*stensv1(:,kr2,d1)*h0022_1;              
            stenv22h2000(:,d1)=stenv22h2000(:,d1)+wtk.*stensv22(:,kr2,d1)*h2000_1;              
            
            stens4v1v1v1(:,d1) = stens4v1v1v1(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v1_1;
            stens4v1v1v11(:,d1) = stens4v1v1v11(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v11_1;
            stens4v1v1v2(:,d1) = stens4v1v1v2(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v2_1;
            stens4v1v1v22(:,d1) = stens4v1v1v22(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v22_1;
            stens4v1v11v2(:,d1) = stens4v1v11v2(:,d1)+(wtk.*stens4v1v11(:,kr2,d1))*v2_1;            
            stens4v1v11v11(:,d1) = stens4v1v11v11(:,d1)+(wtk.*stens4v1v11(:,kr2,d1))*v11_1;
            stens4v1v1v11(:,d1) = stens4v1v1v11(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v11_1;
            stens4v2v2v2(:,d1) = stens4v2v2v2(:,d1)+(wtk.*stens4v2v2(:,kr2,d1))*v2_1;
            stens4v2v22v22(:,d1) = stens4v2v22v22(:,d1)+(wtk.*stens4v2v22(:,kr2,d1))*v22_1;
            stens4v2v2v22(:,d1) = stens4v2v2v22(:,d1)+(wtk.*stens4v2v2(:,kr2,d1))*v22_1;
            stens4v1v22v22(:,d1) = stens4v1v22v22(:,d1)+(wtk.*stens4v1v22(:,kr2,d1))*v22_1;
            stens4v1v2v2(:,d1) = stens4v1v2v2(:,d1)+(wtk.*stens4v1v2(:,kr2,d1))*v2_1;
            stens4v1v2v22(:,d1) = stens4v1v2v22(:,d1)+(wtk.*stens4v1v2(:,kr2,d1))*v22_1;
            stens4v11v11v2(:,d1) = stens4v11v11v2(:,d1)+(wtk.*stens4v11v11(:,kr2,d1))*v2_1;
            stens4v1v11v22(:,d1) = stens4v1v11v22(:,d1)+(wtk.*stens4v1v11(:,kr2,d1))*v22_1;            
            stens4v11v2v22(:,d1) = stens4v11v2v22(:,d1)+(wtk.*stens4v11v2(:,kr2,d1))*v22_1;
            stens4v11v2v2(:,d1) = stens4v11v2v2(:,d1)+(wtk.*stens4v11v2(:,kr2,d1))*v2_1;
                  
            stens5v1v1v1v11(:,d1) = stens5v1v1v1v11(:,d1)+wtk.*stens5v1v1v1(:,kr2,d1)*v11_1;
            stens5v1v1v11v11(:,d1) = stens5v1v1v11v11(:,d1)+wtk.*stens5v1v1v11(:,kr2,d1)*v11_1;
            stens5v2v2v2v22(:,d1) = stens5v2v2v2v22(:,d1)+wtk.*stens5v2v2v2(:,kr2,d1)*v22_1;            
            stens5v1v2v2v22(:,d1) = stens5v1v2v2v22(:,d1)+wtk.*stens5v1v2v2(:,kr2,d1)*v22_1;
            stens5v1v1v11v2(:,d1) = stens5v1v1v11v2(:,d1)+wtk.*stens5v1v1v11(:,kr2,d1)*v2_1;
            stens5v1v11v2v2(:,d1) = stens5v1v11v2v2(:,d1)+wtk.*stens5v1v11v2(:,kr2,d1)*v2_1;

        end

        fxhessv1h0022(range2)=wtk.*shv1(:,kr2)*h0022_1;
        fxhessv1h2200(range2)=wtk.*shv1(:,kr2)*h2200_1;
        fxhessv1h1210(range2)=wtk.*shv1(:,kr2)*h1210_1;
        fxhessv1h0121(range2)=wtk.*shv1(:,kr2)*h0121_1;        
        fxhessv1h1111(range2)=wtk.*shv1(:,kr2)*h1111_1;        
        fxhessv11h2110(range2)=wtk.*shv11(:,kr2)*h2110_1;        
        fxhessv11h2011(range2)=wtk.*shv11(:,kr2)*h2011_1;
        fxhessv11h1021(range2)=wtk.*shv11(:,kr2)*h1021_1;            
        fxhessv11h3100(range2)=wtk.*shv11(:,kr2)*h3100_1;
        fxhessv2h0022(range2)=wtk.*shv2(:,kr2)*h0022_1;        
        fxhessv2h1012(range2)=wtk.*shv2(:,kr2)*h1012_1;
        fxhessv2h2200(range2)=wtk.*shv2(:,kr2)*h2200_1;        
        fxhessv2h2101(range2)=wtk.*shv2(:,kr2)*h2101_1;
        fxhessv2h1111(range2)=wtk.*shv2(:,kr2)*h1111_1;
        fxhessv22h2110(range2)=wtk.*shv22(:,kr2)*h2110_1;
        fxhessv22h1120(range2)=wtk.*shv22(:,kr2)*h1120_1;
        fxhessv22h1021(range2)=wtk.*shv22(:,kr2)*h1021_1;                
        fxhessv22h0031(range2)=wtk.*shv22(:,kr2)*h0031_1;
        fxhessh0002h1020(range2)=wtk.*shh0002(:,kr2)*h1020_1;                        
        fxhessh0002h0030(range2)=wtk.*shh0002(:,kr2)*h0030_1;
        fxhessh0011h1011(range2)=wtk.*shh0011(:,kr2)*h1011_1;        
        fxhessh0012h1010(range2)=wtk.*shh0012(:,kr2)*h1010_1;
        fxhessh0020h0012(range2)=wtk.*shh0020(:,kr2)*h0012_1;                
        fxhessh0020h1002(range2)=wtk.*shh0020(:,kr2)*h1002_1;
        fxhessh0020h1101(range2)=wtk.*shh0020(:,kr2)*h1101_1;
        fxhessh0021h0011(range2)=wtk.*shh0021(:,kr2)*h0011_1;
        fxhessh0021h1100(range2)=wtk.*shh0021(:,kr2)*h1100_1;
        fxhessh0021h1001(range2)=wtk.*shh0021(:,kr2)*h1001_1;        
        fxhessh0101h1020(range2)=wtk.*shh0101(:,kr2)*h1020_1;
        fxhessh0101h2010(range2)=wtk.*shh0101(:,kr2)*h2010_1;                
        fxhessh0110h1011(range2)=wtk.*shh0110(:,kr2)*h1011_1;                                  
        fxhessh0120h1001(range2)=wtk.*shh0120(:,kr2)*h1001_1;
        fxhessh0200h3000(range2)=wtk.*shh0200(:,kr2)*h3000_1;        
        fxhessh0200h2010(range2)=wtk.*shh0200(:,kr2)*h2010_1;
        fxhessh1001h1110(range2)=wtk.*shh1001(:,kr2)*h1110_1;        
        fxhessh1010h1101(range2)=wtk.*shh1010(:,kr2)*h1101_1;
        fxhessh1010h0111(range2)=wtk.*shh1010(:,kr2)*h0111_1;                                  
        fxhessh1110h0011(range2)=wtk.*shh1110(:,kr2)*h0011_1;
        fxhessh1011h1100(range2)=wtk.*shh1011(:,kr2)*h1100_1;
        fxhessh1100h1110(range2)=wtk.*shh1100(:,kr2)*h1110_1;
        fxhessh1200h1010(range2)=wtk.*shh1200(:,kr2)*h1010_1;
        fxhessh2000h1200(range2)=wtk.*shh2000(:,kr2)*h1200_1;                
        fxhessh2000h0210(range2)=wtk.*shh2000(:,kr2)*h0210_1;
        fxhessh2000h0111(range2)=wtk.*shh2000(:,kr2)*h0111_1;
        fxhessh2001h0110(range2)=wtk.*shh2001(:,kr2)*h0110_1;
        fxhessh2100h0011(range2)=wtk.*shh2100(:,kr2)*h0011_1;
        fxhessh2100h0110(range2)=wtk.*shh2100(:,kr2)*h0110_1;
        fxhessh2100h1100(range2)=wtk.*shh2100(:,kr2)*h1100_1;
  
        
        fxtensv1v1h1200(range2)=wtk.*stenv1v1(:,kr2)*h1200_1;
        fxtensv1v1h0210(range2)=wtk.*stenv1v1(:,kr2)*h0210_1;        
        fxtensv1v1h0111(range2)=wtk.*stenv1v1(:,kr2)*h0111_1;
        fxtensv1v11h2100(range2)=wtk.*stenv1v11(:,kr2)*h2100_1;
        fxtensv1v11h1011(range2)=wtk.*stenv1v11(:,kr2)*h1011_1;       
        fxtensv1v11h1110(range2)=wtk.*stenv1v11(:,kr2)*h1110_1;                
        fxtensv1v11h0021(range2)=wtk.*stenv1v11(:,kr2)*h0021_1;        
        fxtensv1v2h1200(range2)=wtk.*stenv1v2(:,kr2)*h1200_1;              
        fxtensv1v2h0111(range2)=wtk.*stenv1v2(:,kr2)*h0111_1;                
        fxtensv1v2h1101(range2)=wtk.*stenv1v2(:,kr2)*h1101_1;                
        fxtensv1v2h0012(range2)=wtk.*stenv1v2(:,kr2)*h0012_1;                
        fxtensv1v22h0021(range2)=wtk.*stenv1v22(:,kr2)*h0021_1;                
        fxtensv1v22h0120(range2)=wtk.*stenv1v22(:,kr2)*h0120_1;                
        fxtensv1v22h1110(range2)=wtk.*stenv1v22(:,kr2)*h1110_1;        
        fxtensv11v11h2010(range2)=wtk.*stenv11v11(:,kr2)*h2010_1;        
        fxtensv11v11h3000(range2)=wtk.*stenv11v11(:,kr2)*h3000_1;
        fxtensv11v2h2100(range2)=wtk.*stenv11v2(:,kr2)*h2100_1;        
        fxtensv11v2h2001(range2)=wtk.*stenv11v2(:,kr2)*h2001_1;        
        fxtensv11v2h1011(range2)=wtk.*stenv11v2(:,kr2)*h1011_1;        
        fxtensv11v22h1020(range2)=wtk.*stenv11v22(:,kr2)*h1020_1;        
        fxtensv11v22h2010(range2)=wtk.*stenv11v22(:,kr2)*h2010_1;        
        fxtensv2v2h0012(range2)=wtk.*stenv2v2(:,kr2)*h0012_1;
        fxtensv2v2h1002(range2)=wtk.*stenv2v2(:,kr2)*h1002_1;        
        fxtensv2v2h1101(range2)=wtk.*stenv2v2(:,kr2)*h1101_1;
        fxtensv2v22h2100(range2)=wtk.*stenv2v22(:,kr2)*h2100_1;        
        fxtensv2v22h1110(range2)=wtk.*stenv2v22(:,kr2)*h1110_1;        
        fxtensv2v22h1011(range2)=wtk.*stenv2v22(:,kr2)*h1011_1;        
        fxtensv2v22h0021(range2)=wtk.*stenv2v22(:,kr2)*h0021_1;
        fxtensv22v22h0030(range2)=wtk.*stenv22v22(:,kr2)*h0030_1;
        fxtensv22v22h1020(range2)=wtk.*stenv22v22(:,kr2)*h1020_1;              
                
        fxtensv1h1100h1100(range2)=wtk.*stenv1h1100(:,kr2)*h1100_1;
        fxtensv1h0200h2000(range2)=wtk.*stenv1h0200(:,kr2)*h2000_1;        
        fxtensv1h0011h0011(range2)=wtk.*stenv1h0011(:,kr2)*h0011_1;                
        fxtensv1h0020h0002(range2)=wtk.*stenv1h0020(:,kr2)*h0002_1;                               
        fxtensv1h0110h1100(range2)=wtk.*stenv1h0110(:,kr2)*h1100_1;        
        fxtensv1h1010h0200(range2)=wtk.*stenv1h1010(:,kr2)*h0200_1;        
        fxtensv1h0020h0101(range2)=wtk.*stenv1h0020(:,kr2)*h0101_1;                
        fxtensv1h1100h0011(range2)=wtk.*stenv1h1100(:,kr2)*h0011_1;        
        fxtensv1h1001h0110(range2)=wtk.*stenv1h1001(:,kr2)*h0110_1;        
        fxtensv1h0110h0011(range2)=wtk.*stenv1h0110(:,kr2)*h0011_1;                        
        fxtensv1h1010h0101(range2)=wtk.*stenv1h1010(:,kr2)*h0101_1;        
        fxtensv11h0020h1001(range2)=wtk.*stenv11h0020(:,kr2)*h1001_1;                
        fxtensv11h1010h0011(range2)=wtk.*stenv11h1010(:,kr2)*h0011_1;        
        fxtensv11h1001h1010(range2)=wtk.*stenv11h1001(:,kr2)*h1010_1;        
        fxtensv11h2000h0011(range2)=wtk.*stenv11h2000(:,kr2)*h0011_1;           
        fxtensv11h2000h1100(range2)=wtk.*stenv11h2000(:,kr2)*h1100_1;
        fxtensv11h0110h2000(range2)=wtk.*stenv11h0110(:,kr2)*h2000_1;        
        fxtensv11h1010h1100(range2)=wtk.*stenv11h1010(:,kr2)*h1100_1; 
        fxtensv2h0110h1001(range2)=wtk.*stenv2h0110(:,kr2)*h1001_1;        
        fxtensv2h2000h0101(range2)=wtk.*stenv2h2000(:,kr2)*h0101_1;        
        fxtensv2h0011h1100(range2)=wtk.*stenv2h0011(:,kr2)*h1100_1;        
        fxtensv2h1001h1100(range2)=wtk.*stenv2h1001(:,kr2)*h1100_1;                
        fxtensv2h1010h0101(range2)=wtk.*stenv2h1010(:,kr2)*h0101_1;        
        fxtensv2h1100h1100(range2)=wtk.*stenv2h1100(:,kr2)*h1100_1;        
        fxtensv2h0002h0020(range2)=wtk.*stenv2h0002(:,kr2)*h0020_1;                                
        fxtensv2h1001h0011(range2)=wtk.*stenv2h1001(:,kr2)*h0011_1;
        fxtensv2h1010h0002(range2)=wtk.*stenv2h1010(:,kr2)*h0002_1;
        fxtensv2h2000h0200(range2)=wtk.*stenv2h2000(:,kr2)*h0200_1;        
        fxtensv2h0011h0011(range2)=wtk.*stenv2h0011(:,kr2)*h0011_1;
        fxtensv22h0020h0011(range2)=wtk.*stenv22h0020(:,kr2)*h0011_1;        
        fxtensv22h1001h0020(range2)=wtk.*stenv22h1001(:,kr2)*h0020_1;       
        fxtensv22h1010h0011(range2)=wtk.*stenv22h1010(:,kr2)*h0011_1;        
        fxtensv22h1010h1100(range2)=wtk.*stenv22h1010(:,kr2)*h1100_1;              
        fxtensv22h0110h1010(range2)=wtk.*stenv22h0110(:,kr2)*h1010_1;        
        fxtensv22h0020h1100(range2)=wtk.*stenv22h0020(:,kr2)*h1100_1;
        fxtensv22h2000h0110(range2)=wtk.*stenv22h2000(:,kr2)*h0110_1;
                                     
        fxtens4v1v1v1h0200(range2)=wtk.*stens4v1v1v1(:,kr2)*h0200_1;
        fxtens4v1v1v11h0011(range2)=wtk.*stens4v1v1v11(:,kr2)*h0011_1;        
        fxtens4v1v1v11h1100(range2)=wtk.*stens4v1v1v11(:,kr2)*h1100_1;        
        fxtens4v1v1v11h0110(range2)=wtk.*stens4v1v1v11(:,kr2)*h0110_1;        
        fxtens4v1v1v2h0101(range2)=wtk.*stens4v1v1v2(:,kr2)*h0101_1;        
        fxtens4v1v1v2h0200(range2)=wtk.*stens4v1v1v2(:,kr2)*h0200_1;            
        fxtens4v1v1v22h0110(range2)=wtk.*stens4v1v1v22(:,kr2)*h0110_1;        
        fxtens4v1v11v11h1010(range2)=wtk.*stens4v1v11v11(:,kr2)*h1010_1;        
        fxtens4v1v11v11h2000(range2)=wtk.*stens4v1v11v11(:,kr2)*h2000_1;
        fxtens4v1v11v2h0011(range2)=wtk.*stens4v1v11v2(:,kr2)*h0011_1;        
        fxtens4v1v11v2h1001(range2)=wtk.*stens4v1v11v2(:,kr2)*h1001_1;                
        fxtens4v1v11v2h1100(range2)=wtk.*stens4v1v11v2(:,kr2)*h1100_1;        
        fxtens4v1v11v22h0020(range2)=wtk.*stens4v1v11v22(:,kr2)*h0020_1;        
        fxtens4v1v11v22h1010(range2)=wtk.*stens4v1v11v22(:,kr2)*h1010_1;        
        fxtens4v1v2v2h0101(range2)=wtk.*stens4v1v2v2(:,kr2)*h0101_1;        
        fxtens4v1v2v2h0002(range2)=wtk.*stens4v1v2v2(:,kr2)*h0002_1;                
        fxtens4v1v2v22h0011(range2)=wtk.*stens4v1v2v22(:,kr2)*h0011_1;        
        fxtens4v1v2v22h0110(range2)=wtk.*stens4v1v2v22(:,kr2)*h0110_1;          
        fxtens4v1v2v22h1100(range2)=wtk.*stens4v1v2v22(:,kr2)*h1100_1;                
        fxtens4v1v22v22h0020(range2)=wtk.*stens4v1v22v22(:,kr2)*h0020_1;                
        fxtens4v11v11v2h2000(range2)=wtk.*stens4v11v11v2(:,kr2)*h2000_1;        
        fxtens4v11v2v2h1001(range2)=wtk.*stens4v11v2v2(:,kr2)*h1001_1;        
        fxtens4v11v2v22h2000(range2)=wtk.*stens4v11v2v22(:,kr2)*h2000_1;        
        fxtens4v11v2v22h1010(range2)=wtk.*stens4v11v2v22(:,kr2)*h1010_1;        
        fxtens4v2v2v2h0002(range2)=wtk.*stens4v2v2v2(:,kr2)*h0002_1;
        fxtens4v2v2v22h1001(range2)=wtk.*stens4v2v2v22(:,kr2)*h1001_1;               
        fxtens4v2v2v22h1100(range2)=wtk.*stens4v2v2v22(:,kr2)*h1100_1;        
        fxtens4v2v2v22h0011(range2)=wtk.*stens4v2v2v22(:,kr2)*h0011_1;                
        fxtens4v2v22v22h0020(range2)=wtk.*stens4v2v22v22(:,kr2)*h0020_1;        
        fxtens4v2v22v22h1010(range2)=wtk.*stens4v2v22v22(:,kr2)*h1010_1;        
                  
        fxtens5v1v1v1v11v11(range2)=wtk.*stens5v1v1v1v11(:,kr2)*v11_1;        
        fxtens5v1v1v11v11v2(range2)=wtk.*stens5v1v1v11v11(:,kr2)*v2_1;        
        fxtens5v1v1v11v2v22(range2)=wtk.*stens5v1v1v11v2(:,kr2)*v22_1;                
        fxtens5v1v11v2v2v22(range2)=wtk.*stens5v1v11v2v2(:,kr2)*v22_1;        
        fxtens5v1v2v2v22v22(range2)=wtk.*stens5v1v2v2v22(:,kr2)*v22_1;        
        fxtens5v2v2v2v22v22(range2)=wtk.*stens5v2v2v2v22(:,kr2)*v22_1;        
       
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end

%computation of a3200
diffh2100 =  T*fxjac*h2100-i*theta1*h2100C+diffh2100.';
a3200 = 1/(12*T^2)*w1'*(fxtens5v1v1v1v11v11.'+fxtens4v1v1v1h0200.'+3*fxtens4v1v11v11h2000.'...
    +6*fxtens4v1v1v11h1100.'+6*fxtensv1h1100h1100.'+3*fxtensv1v1h1200.'+fxtensv11v11h3000.'...
    +6*fxtensv1v11h2100.'+6*fxtensv11h2000h1100.'+3*fxtensv1h0200h2000.'+fxhessh0200h3000.'+2*fxhessv11h3100.'...
    +3*fxhessv1h2200.'+6*fxhessh2100h1100.'+3*fxhessh2000h1200.'-6*alpha11001*diffh2100-12*alpha22001*fxjac*v1)+i*theta1/T*alpha22001/T^2;

%computation of b0032
diffh0021 =  T*fxjac*h0021-i*theta2*h0021C+diffh0021.';
b0032 = 1/(12*T^2)*w2'*(fxtens5v2v2v2v22v22.'+fxtens4v2v2v2h0002.'+3*fxtens4v2v22v22h0020.'...
    +6*fxtens4v2v2v22h0011.'+6*fxtensv2h0011h0011.'+3*fxtensv2v2h0012.'+fxtensv22v22h0030.'...
    +6*fxtensv2v22h0021.'+6*fxtensv22h0020h0011.'+3*fxtensv2h0002h0020.'+fxhessh0002h0030.'+2*fxhessv22h0031.'...
    +3*fxhessv2h0022.'+6*fxhessh0021h0011.'+3*fxhessh0020h0012.'-6*alpha00111*diffh0021-12*alpha00221*fxjac*v2)+i*theta2/T*alpha00221/T^2;

%computation of a1022
diffh1011 =  T*fxjac*h1011-i*theta1*h1011C+diffh1011.';
a1022 = 1/(4*T^2)*w1'*(fxtens5v1v2v2v22v22.'+fxtens4v1v22v22h0020.'+fxtens4v1v2v2h0002.'+2*fxtens4v2v2v22h1001.'...
    +2*fxtens4v2v22v22h1010.'+4*fxtens4v1v2v22h0011.'+2*fxtensv1v22h0021.'+fxtensv1h0020h0002.'+2*fxtensv1v2h0012.'...
    +2*fxtensv22h1001h0020.'+4*fxtensv2h1001h0011.'+2*fxtensv2h1010h0002.'+fxtensv2v2h1002.'+4*fxtensv22h1010h0011.'...
    +4*fxtensv2v22h1011.'+fxtensv22v22h1020.'+2*fxtensv1h0011h0011.'+fxhessv1h0022.'+2*fxhessh0021h1001.'+fxhessh0020h1002.'...
    +2*fxhessh0012h1010.'+4*fxhessh0011h1011.'+2*fxhessv2h1012.'+fxhessh0002h1020.'+2*fxhessv22h1021.'...
    -4*alpha00111*diffh1011-4*alpha00221*fxjac*v1)+i*theta1/T*alpha00221/T^2;

%computation of b2210
diffh1110 =  T*fxjac*h1110-i*theta2*h1110C+diffh1110.';
b2210 = 1/(4*T^2)*w2'*(fxtens5v1v1v11v11v2.'+fxtens4v11v11v2h2000.'+fxtens4v1v1v2h0200.'+2*fxtens4v1v1v11h0110.'...
    +2*fxtens4v1v11v11h1010.'+4*fxtens4v1v11v2h1100.'+2*fxtensv11v2h2100.'+fxtensv2h2000h0200.'+2*fxtensv1v2h1200.'...
    +2*fxtensv11h0110h2000.'+4*fxtensv1h0110h1100.'+2*fxtensv1h1010h0200.'+fxtensv1v1h0210.'+4*fxtensv11h1010h1100.'...
    +4*fxtensv1v11h1110.'+fxtensv11v11h2010.'+2*fxtensv2h1100h1100.'+fxhessv2h2200.'+2*fxhessh2100h0110.'+fxhessh2000h0210.'...
    +2*fxhessh1200h1010.'+4*fxhessh1100h1110.'+2*fxhessv1h1210.'+fxhessh0200h2010.'+2*fxhessv11h2110.'-4*alpha11001*diffh1110...
    -4*alpha22001*fxjac*v2)+i*theta2/T*alpha22001/T^2;
b22101 = T^2*b2210;

%computation of a2111
a2111 = 1/(2*T^2)*w1'*(fxtens5v1v1v11v2v22.'+fxtens4v1v1v2h0101.'+fxtens4v1v1v11h0011.'+fxtens4v1v1v22h0110.'...
    +2*fxtens4v1v11v22h1010.'+2*fxtens4v1v2v22h1100.'+fxtens4v11v2v22h2000.'+2*fxtens4v1v11v2h1001.'+2*fxtensv1h1001h0110.'...
    +fxtensv2v22h2100.'+2*fxtensv2h1001h1100.'+2*fxtensv1h1100h0011.'+2*fxtensv1h1010h0101.'+fxtensv11v22h2010.'+fxtensv11v2h2001.'...
    +2*fxtensv1v22h1110.'+2*fxtensv22h1010h1100.'+2*fxtensv1v2h1101.'+fxtensv2h2000h0101.'+2*fxtensv11h1001h1010.'...
    +fxtensv11h2000h0011.'+fxtensv22h2000h0110.'+2*fxtensv1v11h1011.'+fxtensv1v1h0111.'+fxhessv2h2101.'+fxhessh2100h0011.'...
    +2*fxhessv1h1111.'+fxhessv11h2011.'+fxhessh0101h2010.'+fxhessh2001h0110.'+fxhessh2000h0111.'+fxhessv22h2110.'...
    +2*fxhessh1011h1100.'+2*fxhessh1010h1101.'+2*fxhessh1001h1110.'-2*alpha11111*fxjac*v1-2*alpha11001*diffh1011...
    -alpha00111*diffh2100)+i*theta1/T*alpha11111/T^2;
a21111 = T^2*a2111;

%computation of b1121
b1121 = 1/(2*T^2)*w2'*(fxtens5v1v11v2v2v22.'+fxtens4v1v2v2h0101.'+fxtens4v2v2v22h1100.'+fxtens4v11v2v2h1001.'...
    +2*fxtens4v11v2v22h1010.'+2*fxtens4v1v11v2h0011.'+fxtens4v1v11v22h0020.'+2*fxtens4v1v2v22h0110.'+2*fxtensv2h0110h1001.'...
    +fxtensv1v11h0021.'+2*fxtensv1h0110h0011.'+2*fxtensv2h0011h1100.'+2*fxtensv2h1010h0101.'+fxtensv11v22h1020.'...
    +fxtensv1v22h0120.'+2*fxtensv11v2h1011.'+2*fxtensv11h1010h0011.'+2*fxtensv1v2h0111.'+fxtensv1h0020h0101.'...
    +2*fxtensv22h0110h1010.'+fxtensv22h0020h1100.'+fxtensv11h0020h1001.'+2*fxtensv2v22h1110.'+fxtensv2v2h1101.'...
    +fxhessv1h0121.'+fxhessh0021h1100.'+2*fxhessv2h1111.'+fxhessv22h1120.'+fxhessh0101h1020.'+fxhessh0120h1001.'...
    +fxhessh0020h1101.'+fxhessv11h1021.'+2*fxhessh1110h0011.'+2*fxhessh1010h0111.'+2*fxhessh0110h1011.'-2*alpha11111*fxjac*v2...
    -2*alpha00111*diffh1110-alpha11001*diffh0021)+i*theta2/T*alpha11111/T^2;
b11211 = T^2*b1121;
end


%transformation formulas
G2100 = a2100-i*theta1/T*alpha1100;
G1011 = a1011-i*theta1/T*alpha0011;
H1110 = b1110-i*theta2/T*alpha1100;
H0021 = b0021-i*theta2/T*alpha0011;
if higherorderderivatives == 1
G3200 = a3200-i*theta1/T*alpha2200+(-a2100+i*theta1/T*alpha1100)*alpha1100;
G2111 = -i*theta1/T*alpha1111+a2111+(-a1011+i*theta1/T*alpha0011)*alpha1100+(-a2100+i*theta1/T*alpha1100)*alpha0011;
G1022 = a1022+(-a1011+i*theta1/T*alpha0011)*alpha0011-i*theta1/T*alpha0022;
H2210 = b2210+(-b1110+i*theta2/T*alpha1100)*alpha1100-i*theta2/T*alpha2200;
H1121 = -i*theta2/T*alpha1111+b1121+(-b1110+i*theta2/T*alpha1100)*alpha0011+(-b0021+i*theta2/T*alpha0011)*alpha1100;
H0032 = b0032-i*theta2/T*alpha0022+(-b0021+i*theta2/T*alpha0011)*alpha0011;
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

% -------------------------------------------------------------
function [x,p,T] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global lds

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+1:lds.PeriodIdx+2);
x = x0(lds.coords);
T = x0(lds.PeriodIdx);


% -------------------------------------------------------------

function [jac1,jac2,jac3] = BVP_jac(BVP_jac_f,BVP_jac_bc1,x,p,T,theta1,theta2,pars,nc)
global lds 

ups = reshape(x,lds.nphase,lds.tps);
p = num2cell(p);

jac1 = spalloc(lds.ncoords,lds.ncoords,(lds.ncol+4)*lds.nphase);

range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
jac2 = jac1;
eitheta1 = 1i*theta1*lds.pwwt;
eitheta2 = 1i*theta2*lds.pwwt;
for jk = lds.tsts
  % value of polynomial on each collocation point
  xp = ups(:,range0)*lds.wt;

  % evaluate part of Jacobian matrix
  jac3(range1,range2) = feval(BVP_jac_f,lds.func,xp,p,T,jk);
  jac1(range1,range2) = jac3(range1,range2)+eitheta1;
  jac2(range1,range2) = jac3(range1,range2)+eitheta2;

  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end

% boundary conditions
range = (lds.tps-1)*lds.nphase + (lds.phases);
jac1(range,[lds.phases range]) = feval(BVP_jac_bc1);
jac2(range,[lds.phases range]) = jac1(range,[lds.phases range]);
jac3(range,[lds.phases range]) = jac1(range,[lds.phases range]);

