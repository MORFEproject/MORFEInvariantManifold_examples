function result=nf_LPNS(x)
%
% calculates LPNS normal form coefficient
%
global lds

higherorderderivatives = 0;

alpha200 = 0;
alpha011 = 0;

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

%computation f1
% function
range1 = lds.cols_p1;
range2 = lds.phases;
for jk=lds.tsts
    xpp = ups(:,range1)*lds.wt;
    for c=lds.cols                    
        xt = xpp(:,c);
        f1(range2) = feval(lds.func, 0,xt, pt{:});
        range2 = range2 + lds.nphase;
    end   
    range1 = range1 + lds.ncol;
end
f1 = f1';

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

%computation of v1
b = []; b(lds.ncoords+1)=1; b=b';
jace=[Jjac2,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

Jjac2=[Jjac2 p;ic 0];
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*f1; rl(lds.ncoords+1) = 0;
v1 = Jjac2\rl';
v1 = v1(lds.coords);
Jjac2(lds.ncoords+1,lds.coords) = q';

%computation of phi*
phi = Jjac2'\b;
phi = phi(1:lds.ncoords-lds.nphase)';

%rescale phi
v1ps = reshape(v1,lds.nphase,lds.tps);
range1 = lds.cols_p1;
range2 = 1 : lds.ncol;
for j=lds.tsts
    v1Cps(:,range2) = v1ps(:,range1)*lds.wt;
    range1 = range1 + lds.ncol;
    range2 = range2 + lds.ncol;
end
v1C = reshape(v1Cps,lds.ncoords-lds.nphase,1);
phi = phi/(phi*v1C);

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

%computation of v2* (up to scaling)
wext = LJ'\(UJ'\b);
w2 = wext(1:lds.ncoords-lds.nphase);

phips = reshape(phi,lds.nphase,lds.tps-1);
v2ps = reshape(v2,lds.nphase,lds.tps);
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
range3 = lds.cols_p1;
wt = lds.wt';
for j=lds.tsts
    v2Cps(:,range2) = v2ps(:,range3)*lds.wt;
    p = phips(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
    range3 = range3 + lds.ncol;
end
v2C = reshape(v2Cps,lds.ncoords-lds.nphase,1);

%computation of w1
p = Jjac2(lds.coords,lds.ncoords+1);
Jjac2(:,lds.ncoords+1) = [v1C;zeros(lds.nphase+1,1)];
rl(1:lds.ncoords) = T*ic; 
rl(lds.ncoords+1) = 0;
w1 = Jjac2'\rl';
w1 = w1(1:lds.ncoords-lds.nphase);
Jjac2(lds.coords,lds.ncoords+1) = p;

%rescale w2
ic = w2'*v2C;
w2 = w2/ic';

%computation of a200
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
w22ps = conj(w2ps);
icw1 = zeros(1,lds.ncoords);
icw22 = zeros(1,lds.ncoords);
wt = lds.wt';
for j=lds.tsts
    % value of polynomial on each collocation point
    xpp(:,range6) = ups(:,range1)*lds.wt;
    xp = xpp(:,range6);        
    v1_1 = v1(range4);
    v2_1 = v2(range4);
    v22_1 = conj(v2_1);
    p = w1ps(:,range6)*wt;
    icw1(range4) = icw1(range4)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
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

a200 = 1/2*phi*(fxhessv1v1'+2*fxjac*v1);

%computation of h200
Jjac2(lds.ncoords+1,lds.coords) = icw1;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v1'+2*T*fxjac*v1-2*T*a200*v1C-2*T*alpha200*f1+2*T*f1; rl(lds.ncoords+1) = 0;
diffh200 = rl(1:lds.ncoords-lds.nphase);
h200 = Jjac2\rl';
h200 = h200(lds.coords);

%computation of h020
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
range3 = lds.cols;
twitheta = 2*1i*theta*lds.pwwt;
J = Jjac;
for jk  = lds.tsts
  % evaluate part of Jacobian matrix
  J(range1,range2) = bordBVP_R2_f(lds.func,xpp(:,range3),pt,T,jk)+twitheta;  
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
  range3 = range3 + lds.ncol;
end
% remove borders
J(:,end)=[];J(end,:)=[];
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv2v2; rl(lds.ncoords) = 0;
h020 = J\rl.';

%computation of b110
b110 = w2'*(fxhessv1v2.'+fxjac*v2)-i*theta/T;

%computation of h110
Jjac(lds.ncoords+1,lds.coords) = icw22;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v2.'+T*fxjac*v2-T*b110*v2C-i*theta*v2C; rl(lds.ncoords+1) = 0;
diffh110 = rl(1:lds.ncoords-lds.nphase) ;
h110 = Jjac\rl.';
h110 = h110(lds.coords);

%computation of a0111
a0111 = phi*real(fxhessv2v22)';
a011 = a0111/T;

%computation of h011
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxhessv2v22)'-T*a0111*v1C-T^2*alpha011*f1; rl(lds.ncoords+1) = 0;
diffh011 = rl(1:lds.ncoords-lds.nphase);
h011 = Jjac2\rl.';
h011 = h011(lds.coords);

if higherorderderivatives == 1
% function
range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
h200ps = reshape(h200,lds.nphase,lds.tps);
h011ps = reshape(h011,lds.nphase,lds.tps);
h110ps = reshape(h110,lds.nphase,lds.tps);
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    h200Cps(:,range6) = h200ps(:,range1)*lds.wt;
    h011Cps(:,range6) = h011ps(:,range1)*lds.wt;
    h110Cps(:,range6) = h110ps(:,range1)*lds.wt;
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
        fxhessv2h101(range2)=wtk.*shv2(:,kr2)*h101_1;
        fxhessv1h011(range2)=wtk.*shv1(:,kr2)*h011_1;
        fxhessv2h011(range2)=wtk.*shv2(:,kr2)*h011_1;
        fxhessv22h020(range2)=wtk.*shv22(:,kr2)*h020_1;        
        fxtensv1v1v1(range2)=wtk.*stenv1v1(:,kr2)*v1_1;
        fxtensv1v1v2(range2)=wtk.*stenv1v1(:,kr2)*v2_1;
        fxtensv1v2v22(range2)=wtk.*stenv1v2(:,kr2)*v22_1;
        fxtensv2v2v22(range2)=wtk.*stenv2v2(:,kr2)*v22_1;
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end
h200C = reshape(h200Cps,lds.ncoords-lds.nphase,1);
h011C = reshape(h011Cps,lds.ncoords-lds.nphase,1);
h110C = reshape(h110Cps,lds.ncoords-lds.nphase,1);

%computation of a300
diffh200 = T*fxjac*h200+diffh200.';
a300 = 1/6*phi*(fxtensv1v1v1.'+3*fxhessv1h200.'+3*diffh200-6*a200*h200C-6*alpha200*fxjac*v1);

%computation of b210
diffh110 = T*fxjac*h110-i*theta*h110C+diffh110.';
b210 = 1/2*w2'*(fxtensv1v1v2.'+fxhessv2h200.'+2*fxhessv1h110.'+2*diffh110-2*alpha200*fxjac*v2)+i*theta/T*alpha200;

%computation of b021
b021 = 1/(2*T)*w2'*(fxtensv2v2v22.'+fxhessv22h020.'+2*fxhessv2h011.'-2*T*alpha011*fxjac*v2)+i*theta/T*alpha011;

%computation of a111
diffh011 = T*fxjac*h011+diffh011.';
a111 = 1/T*phi*(real(fxtensv1v2v22)'+2*real(fxhessv2h101)'+fxhessv1h011'+diffh011-alpha011*fxjac*v1*T-2*real(b110)*h011C-a0111*h200C);
end 

G200 = a200;
G011 = a011;
H110 = b110+i*theta/T;
if higherorderderivatives == 1
G111 = a111+a011;
G300 = a300+a200;
H210 = b210-i*theta/T*alpha200+b110+i*theta/T;
H021 = b021-i*theta/T*alpha011;
end 

s = sign(G200*G011);
theta1 = real(H110)/G200;
if higherorderderivatives == 1
E = real(H210+H110*(real(H021)/G011-3*G300/(2*G200)+G111/(2*G011))-H021*G200/G011);
else
E = NaN;
end 

result = [s theta1 E];

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

function [jac1,jac2] = BVP_jac(BVP_jac_f,BVP_jac_bc,x,p,T,theta,pars,nc)
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
jac1(range,[lds.phases range]) = feval(BVP_jac_bc);
jac2(range,[lds.phases range]) = jac1(range,[lds.phases range]);