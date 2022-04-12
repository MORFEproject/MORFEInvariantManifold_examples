function result=nf_R1(x)
%
% calculates R1 normal form coefficient
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
jac(range,[lds.phases range1]) = bordBVP_LPC_bc1;

jac(lds.ncoords+1,lds.coords) = bordBVP_CPC_bc2(lds.func,ups,pt);


J = [jac(lds.coords,lds.coords) rand(lds.ncoords,1);rand(1,lds.ncoords) 0];
b=zeros(lds.ncoords+1,1);
b(lds.ncoords+1,1)=1;
p=J'\b;
p=p(1:lds.ncoords,1);p1=p/norm(p);
q=J\b;
q=q(1:lds.ncoords,1);q1=q/norm(q);
J = [jac(:,lds.coords) [p1;0]];
rl = -T*jac(1:(lds.ncoords-lds.nphase),lds.ncoords+1);
rl(lds.ncoords+1) = 0;
vext = J\rl;
v1 = vext(lds.coords);

%compute v21M
%convert v11M to v11C
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

rl = -T*v1C;
rl(lds.ncoords+1) = 0;

v2 = J\rl;
v2 = v2(lds.coords);

%computation of phi
J = [jac(lds.coords,lds.coords) p1;q1' 0];
phi = J'\b;
phi = phi(1:lds.ncoords-lds.nphase);

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
ic = ic*v2;
phi = phi/ic;

% function
range1 = lds.cols_p1;
range2 = lds.phases;
range3 = lds.col_coords;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
v2ps = reshape(v2,lds.nphase,lds.tps);
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    v2Cps(:,range6) = v2ps(:,range1)*lds.wt;    
    v1_1 = v1(range4);
    v2_1 = v2(range4);    
    range5 = lds.phases;             
    % evaluate function value on each collocation point
    for c=lds.cols               
        xt = xp(:,c);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        sjac  = cjac(lds.func,lds.Jacobian,xt,pt,lds.ActiveParams);
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);                    
        sysjac(range5,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',sjac);
        for d=lds.phases
            shv1(:,d) = (wtk.*hess(:,kr2,d))*v1_1;
        end
        fxhessv1v1(range2)=wtk.*shv1(:,kr2)*v1_1;
        fxhessv1v2(range2)=wtk.*shv1(:,kr2)*v2_1;
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end    
    fxjac(range3,range4)=sysjac;
    range1 = range1 + lds.ncol;
    range3 = range3 + lds.ncol_coord;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end
v2C = reshape(v2Cps,lds.ncoords-lds.nphase,1);

%computation of a1
a1 = 1/2*phi'*(2*fxjac*v1+fxhessv1v1');

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

J(1:lds.ncoords-lds.nphase,lds.ncoords+1) = v2C;
J(lds.ncoords-lds.nphase+1:lds.ncoords+1,lds.ncoords+1) = zeros(lds.nphase+1,1);
rl = -T*ic;
rl = [rl 0];
wext = J'\rl';
w1 = wext(1:lds.ncoords-lds.nphase);
b1 = phi'*(fxhessv1v2'+fxjac*v2)+w1'*(2*fxjac*v1+fxhessv1v1');
result=a1*b1;

% ---------------------------------------------------------------
function [x,p,T] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global lds

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+1:lds.PeriodIdx+2);
x = x0(lds.coords);
T = x0(lds.PeriodIdx);


