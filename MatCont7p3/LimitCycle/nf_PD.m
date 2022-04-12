function result = nf_PD(x)
%
% calculates flip normal form coefficient
%
global lds
[xt,p,T] = rearr(x);
pt = num2cell(p);
J = BVP_jac('BVP_PD_jac',xt,p,T,1,1);
[LJ,UJ] = lu(J);
b = []; b(lds.ncoords+1)=1; b=b';
wext = LJ'\(UJ'\b);
vext = UJ\(LJ\b);
% %compute borders
v1 = vext(1:lds.ncoords);
vps = reshape(v1,lds.nphase,lds.tps);
% rescale vext
ficd = dot(vps,vps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v1 = v1/sqrt(ic);
% rescale wext

w1 = wext(1:lds.ncoords-lds.nphase)';
w2 = reshape(w1,lds.nphase,lds.tps-1);
ic2 = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
wt = lds.wt';
for j=lds.tsts
    pw =w2(:,range2)*wt;
    ic2(range1) = ic2(range1)+reshape(pw,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end
%norm v*
ic3 = norm(w1,1);
w1=w1/ic3;

%calculate psi*

% change boundary conditions
range = (lds.tps-1)*lds.nphase + (lds.phases);
J(range,[lds.phases range]) = BVP_LC1_jac_bc;
% remove borders
J(:,end)=[];J(end,:)=[];

%compute borders


jace=[J,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;p=p(1:lds.ncoords,1);p=p/norm(p);

%add borders
J=[J p;q' 0];

psi = J'\b;

ups = reshape(xt,lds.nphase,lds.tps);

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
    v3 = v1(range4);
    range5 = lds.phases;         
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        sten = zeros(lds.nphase,lds.nphase);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        sjac  = cjac(lds.func,lds.Jacobian,xt,pt,lds.ActiveParams);
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        f1(range2) = feval(lds.func, 0,  xt, pt{:});
        sysjac(range5,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',sjac);

        for d1=lds.phases
            sh(:,d1) = (wtk.*hess(:,kr2,d1))*v3;
            for d2=lds.phases
                stens(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v3;
            end    
            sten(:,d1)=sten(:,d1)+wtk.*stens(:,kr2,d1)*v3;
        end
        fxhess(range2)=wtk.*sh(:,kr2)*v3;
        fxtens(range2)=wtk.*sten(:,kr2)*v3;
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    fxjac(range3,range4)=sysjac;
    range1 = range1+lds.ncol;
    range3 = range3 + lds.ncol_coord;
    range4 = range4 + lds.ncol_coord;
end
%rescale psi
psi = psi(1:lds.ncoords-lds.nphase)';
psi = psi*(1/(2*psi*f1'));
%computation of a
a=(psi*fxhess');

%computation of h2
psi1 = reshape(psi,lds.nphase,lds.tps-1);
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
wt = lds.wt';
for j=lds.tsts
    p =psi1(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end
%add integral constraint h2
J(lds.ncoords+1,[lds.coords]) =ic;

b = [];
b(1:lds.ncoords-lds.nphase) = fxhess-2*a*f1; b(lds.ncoords+1) = 0;
h2 = J\b';
range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    v3 = v1(range4);
    h3 = h2(range4);
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        for d1=lds.phases
            sh(:,d1) = (wtk.*hess(:,kr2,d1))*v3;
        end
        fxhessh2(range2) = wtk.*sh(:,kr2)*h3;
        range2 = range2+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
end
result = 1/ic3*ic2*v1*1/3*(w1*(1/T*fxtens'+3*fxhessh2'-6/T*a*fxjac*v1));

% ---------------------------------------------------------------
function [x,p,T] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global lds

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+1);
x = x0(lds.coords);
T = x0(lds.PeriodIdx);

% -------------------------------------------------------------

function jac = BVP_jac(BVP_func,x,p,T,pars,nc)
global lds 
 
p2 = num2cell(p);
jac = feval(BVP_func,lds.func,x,p,T,pars,nc,lds,p2,lds.Jacobian,lds.ActiveParams,lds.JacobianP); 

