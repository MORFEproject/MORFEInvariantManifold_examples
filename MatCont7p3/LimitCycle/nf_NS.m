function result = nf_NS(x0)
%
% calculates NS normal form coefficient
%
global lds 
[xt,p,T] = rearr(x0);pt = num2cell(p);
%computing theta

d = lds.multipliers;
smallest_sum = Inf;
    for jk=1:lds.nphase-1
        [val,idx] = min(abs(d(jk+1:lds.nphase)*d(jk)-1));
        if val < smallest_sum

            idx2 = jk+idx;
            smallest_sum = val;
        end
    end
theta = abs(angle(d(idx2)));

if imag(d(idx2)) == 0 
 % debug('Neutral saddle\n');    
  result='Neutral saddle';
  return;
end

%computing vext and wext

[Jjac,Jjac2] = BVP_jac('bordBVP_PD_jac_f','BVP_LC1_jac_bc',xt,p,T,theta,1,1);

b = []; b(lds.ncoords+1)=1; b=b';

jace=[Jjac,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;p=p(1:lds.ncoords,1);p=p/norm(p);

Jjac=[Jjac p;q.' 0];




[LJ,UJ] = lu(Jjac);

wext = LJ'\(UJ'\b);
vext = UJ \(LJ\b);
ups = reshape(xt,lds.nphase,lds.tps);
%rescale vext
v1 = vext(1:lds.ncoords);
vps = reshape(v1,lds.nphase,lds.tps);
ficd = dot(vps,vps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v1 = v1/sqrt(ic);
% rescale wext
w1 = wext(1:lds.ncoords-lds.nphase).';
ic3 = norm(w1,1);
w1 = w1/ic3;


%calculate psi*

jace=[Jjac2,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;p=p(1:lds.ncoords,1);p=p/norm(p);

Jjac2=[Jjac2 p;q.' 0];



psi = Jjac2'\b;
% function
range1 = lds.cols_p1;
range2 = lds.phases;
range3 = lds.col_coords;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
vps = reshape(v1,lds.nphase,lds.tps);

for jk=lds.tsts
    % value of polynomial on each collocation point
    xp(:,range6) = ups(:,range1)*lds.wt;
    vp(range3) = vps(:,range1)*lds.wt;
    v3 = v1(range4);
    v4 = conj(v3);
    range5 = lds.phases;         
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,range6(c));
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
        fxhess(range2) = wtk.*sh(:,kr2)*v3;
        fxhessc(range2)= wtk.*sh(:,kr2)*v4;
        fxtens(range2) = wtk.*sten(:,kr2)*v4;
        range2 = range2 + lds.nphase;
        range5 = range5 + lds.nphase;
    end   
    fxjac(range3,range4) = sysjac;
    range1 = range1 + lds.ncol;
    range3 = range3 + lds.ncol_coord;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end
alpha = dot(w1,vp);
%rescale psi
psi = psi(1:lds.ncoords-lds.nphase)';
psi = psi/(psi*f1');
%computation of a
a = dot(psi,real(fxhessc));
%calculate h20

range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
range3 = lds.cols;
twitheta = 2*1i*theta*lds.pwwt;
for jk  = lds.tsts
  % evaluate part of Jacobian matrix
  Jjac(range1,range2) = bordBVP_PD_jac_f(lds.func,xp(:,range3),pt,T,jk)+twitheta;  
  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
  range3 = range3 + lds.ncol;
end
% remove borders
Jjac(:,end)=[];Jjac(end,:)=[];
b=[];
b(1:lds.ncoords-lds.nphase) = fxhess; b(lds.ncoords) = 0;
h20 =Jjac\b.';
%computation of h11
psi1 = reshape(psi,lds.nphase,lds.tps-1);
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = lds.cols;
wt = lds.wt';
for jk=lds.tsts
    pw = psi1(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(pw,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end
%add integral constraint h11
Jjac2(lds.ncoords+1,[lds.coords]) =ic;
b = [];
b(1:lds.ncoords-lds.nphase) = real(fxhessc)-a*f1; b(lds.ncoords+1) = 0;
h11 = Jjac2\b.';
range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
for jk=lds.tsts
    % value of polynomial on each collocation point
    v2 = v1(range4);
    v3 = conj(v1(range4));
    h113 = h11(range4);
    h203 = h20(range4);

    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,range6(c));
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams); 
        for d1=lds.phases
            sh11(:,d1) = (wtk.*hess(:,kr2,d1))*h113;
            sh20(:,d1) = (wtk.*hess(:,kr2,d1))*h203;
        end
        fxhessh11(range2) = wtk.*sh11(:,kr2)*v2;
        fxhessh20(range2) = wtk.*sh20(:,kr2)*v3;
        range2 = range2+lds.nphase;
    end   
    range1 = range1 + lds.ncol;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end
result =conj(alpha)*(1/2*dot(w1,(1/T)*(fxtens.')+2*(fxhessh11.')+(fxhessh20.'))-(a/T)*dot(w1,fxjac*v1))+abs(alpha)^2*(1i*a*theta)/(T^2);
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

function [jac1,jac2] = BVP_jac(BVP_jac_f,BVP_jac_bc,x,p,T,theta,pars,nc)
global lds 

ups = reshape(x,lds.nphase,lds.tps);
p = num2cell(p);

jac1 = spalloc(lds.ncoords,lds.ncoords,(lds.ncol+4)*lds.nphase);
%jac = zeros(lds.ncoords+1,lds.ncoords+length(p)-1);

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