function result=nf_R3(x)
%
% calculates R3 normal form coefficient
%
global lds

lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);

[xt,p,T] = rearr(x);


pt = num2cell(p);
pars1 = lds.ncoords+1;
jac = spalloc(lds.ncoords,lds.ncoords,lds.ncol*(lds.ncol+1)*lds.nphase^2*lds.ntst);

ups = reshape(xt,lds.nphase,lds.tps);
%compute v11M
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
jac(range,[lds.phases range1]) = bordBVP_R3_min2pi3_bc1;

b = []; b(lds.ncoords+1)=1; b=b';

jace=[jac(lds.coords,lds.coords),rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

jace = [jac(lds.coords,lds.coords) p;q' 0];
vext = jace\b;
%rescale vext
v1 = vext(1:lds.ncoords);
vps = reshape(v1,lds.nphase,lds.tps);
ficd = dot(vps,vps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v1 = v1/sqrt(ic);

wext = jace'\b;
w1 = wext(1:lds.ncoords-lds.nphase);

v1ps = reshape(v1,lds.nphase,lds.tps);
range0 = lds.cols_p1;
range1 = lds.cols;
for j=lds.tsts
    %value of polynomial on each collocation point
    v1Cps(:,range1) = v1ps(:,range0)*lds.wt;
    range0 = range0 + lds.ncol;
    range1 = range1 + lds.ncol;
end
v1C = reshape(v1Cps,lds.ncoords-lds.nphase,1);

ic = w1'*v1C;
w1 = w1/conj(ic);

% function
range1 = lds.cols_p1;
range2 = lds.phases;
range3 = lds.col_coords;
range4 = lds.cols_p1_coords;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;   
    v1_1 = v1(range4);
    v11_1 = conj(v1_1);    
    range5 = lds.phases;             
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        sjac = cjac(lds.func,lds.Jacobian,xt,pt,lds.ActiveParams);
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);                    
        sysjac(range5,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',sjac);
        for d=lds.phases
            shv1(:,d) = (wtk.*hess(:,kr2,d))*v1_1;
            shv11(:,d) = (wtk.*hess(:,kr2,d))*v11_1;
        end
        fxhessv1v1(range2)=wtk.*shv1(:,kr2)*v1_1;
        fxhessv1v11(range2)=wtk.*shv1(:,kr2)*v11_1;        
        fxhessv11v11(range2)=wtk.*shv11(:,kr2)*v11_1;
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end      
    fxjac(range3,range4)=sysjac;
    range1 = range1 + lds.ncol;
    range3 = range3 + lds.ncol_coord;
    range4 = range4 + lds.ncol_coord;
end

%computation of b1
b1 = 1/2*w1'*fxhessv11v11.';

% computation of h20
w1ps = reshape(w1,lds.nphase,lds.tps-1);
w11ps = conj(w1ps);
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
wt = lds.wt';
for j=lds.tsts
    p = w11ps(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end

jace(lds.ncoords+1,lds.coords) = ic;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv11v11.'-2*b1*T*v1C; rl(lds.ncoords+1) = 0;
h20 = jace\rl.';
h20 = conj(h20(lds.coords));

%computation of phi
% boundary conditions
range  = (lds.tps-1)*lds.nphase+ (lds.phases);
range1 = lds.ncoords-lds.nphase+lds.phases;
jac(range,[lds.phases range1]) = bordBVP_LPC_bc1;
b = []; b(lds.ncoords+1)=1; b=b';

jace=[jac(lds.coords,lds.coords),rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

jace = [jac(lds.coords,lds.coords) p;q' 0];
rl = [];rl(lds.ncoords+1) = 1;
phi = jace'\rl';
phi = phi(1:lds.ncoords-lds.nphase);

ic= -phi'*jac(1:lds.ncoords-lds.nphase,lds.ncoords+1);
phi = phi/ic;

%computation of alpha1
alpha11 = phi'*real(fxhessv1v11)';

%computation of h11
phips = reshape(phi,lds.nphase,lds.tps-1);
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
wt = lds.wt';
for j=lds.tsts
    p = phips(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end

jace(lds.ncoords+1,lds.coords) = ic;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxhessv1v11)'+alpha11*T*jac(1:lds.ncoords-lds.nphase,lds.ncoords+1); rl(lds.ncoords+1) = 0;
h11 = jace\rl';
h11 = h11(lds.coords);

% function
range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
h20ps = reshape(h20,lds.nphase,lds.tps);
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    h20Cps(:,range6) = h20ps(:,range1)*lds.wt;
    v1_1 = v1(range4);
    v11_1 = conj(v1_1);    
    h20_1 = h20(range4);
    h11_1 = h11(range4);
    range5 = lds.phases;             
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        stenv1v1 = zeros(lds.nphase,lds.nphase);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);                    
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        for d=lds.phases
            shv1(:,d) = (wtk.*hess(:,kr2,d))*v1_1;
            shv11(:,d) = (wtk.*hess(:,kr2,d))*v11_1;
            for d2=lds.phases
                stensv1(:,d2,d) = (wtk.*tens(:,kr2,d,d2))*v1_1;
            end    
            stenv1v1(:,d)=stenv1v1(:,d)+wtk.*stensv1(:,kr2,d)*v1_1;
        end
        fxhessv1h11(range2)=wtk.*shv1(:,kr2)*h11_1;
        fxhessv11h20(range2)=wtk.*shv11(:,kr2)*h20_1;        
        fxtensv1v1v11(range2)=wtk.*stenv1v1(:,kr2)*v11_1;
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;        
    end      
    range1 = range1 + lds.ncol;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end
h20C = reshape(h20Cps,lds.ncoords-lds.nphase,1);

c = 1/(2*T)*w1'*(fxtensv1v1v11.'+2*fxhessv1h11.'+fxhessv11h20.'-2*alpha11*fxjac*v1)-1/T*conj(b1)*w1'*conj(h20C);
result = [b1/sqrt(T) real(c)];

% ---------------------------------------------------------------
function [x,p,T] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global lds
p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+1:lds.PeriodIdx+2);
x = x0(lds.coords);
T = x0(lds.PeriodIdx);



