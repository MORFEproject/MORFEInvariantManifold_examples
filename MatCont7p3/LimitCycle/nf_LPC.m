function result=nf_LPC(x)
%
% calculates LPC normal form coefficient
%
global lds

lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);

[xt,p,T] = rearr(x);

pt = num2cell(p);
pars1 = lds.ncoords+1;
jac = spalloc(lds.ncoords,lds.ncoords,(lds.ncol+4)*lds.nphase);

ups = reshape(xt,lds.nphase,lds.tps);
%compute v1M and w1W
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

range = range(lds.nphase)+1;
jac(range,lds.coords) = bordBVP_LPC_bc2(lds.func,ups,pt);

%compute borders
jace=[jac,rand(lds.ncoords+1,1);rand(1,lds.ncoords+1),0];
b=zeros(lds.ncoords+2,1);
b(lds.ncoords+2,1)=1;
q=jace\b;q=q(1:lds.ncoords+1,1);q=q/norm(q);
p=jace'\b;p=p(1:lds.ncoords+1,1);p=p/norm(p);
lds.LPC_psi=p';
lds.LPC_phi=q';

J=[jac lds.LPC_psi';lds.LPC_phi 0];
vext = J\b;
wext = J'\b;
vn = vext(1:lds.ncoords)';
S = vext(lds.ncoords+1);
v1 =(T/S)*vn;
w1 = wext(1:lds.ncoords-lds.nphase)';



ups = reshape(xt,lds.nphase,lds.tps);
vps = reshape(v1',lds.nphase,lds.tps);

% function
range1 = lds.cols_p1;
range2 = lds.phases;
range3 = lds.col_coords;
range4 = lds.cols_p1_coords;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
fxhess = zeros(1,lds.ncoords-lds.nphase);
v2 = zeros(1,lds.ncoords-lds.nphase);
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    vp = vps(:,range1)*lds.wt;
    v3 = v1(range4);
    range5 = lds.phases;         
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        jac  = cjac(lds.func,lds.Jacobian,xt,pt,lds.ActiveParams);
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);                    
        sysjac(range5,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',jac);
        v2(range2) = vp(:,c);
        for d=lds.phases
            sh(:,d) = (wtk.*hess(:,kr2,d))*v3';
        end
        fxhess(range2)=wtk.*sh(:,kr2)*v3';
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    fxjac(range3,range4)=sysjac;
    range1 = range1+lds.ncol;
    range3 = range3 + lds.ncol_coord;
    range4 = range4 + lds.ncol_coord;
end
sysjacv1=fxjac*v1';
result= (w1*(fxhess'+2*sysjacv1))/(2*(w1*v2'));

% ---------------------------------------------------------------
function [x,p,T] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global lds

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+1);
x = x0(lds.coords);
T = x0(lds.PeriodIdx);



