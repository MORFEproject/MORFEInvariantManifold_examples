function result = nf_R2(x)
%
% calculates R2 normal form coefficient
%
global lds

lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);

[xt,p,T] = rearr(x);

pt = num2cell(p);
jac = spalloc(lds.ncoords,lds.ncoords,lds.ncol*(lds.ncol+1)*lds.nphase^2*lds.ntst);
pars1 = lds.ncoords+1;
ups = reshape(xt,lds.nphase,lds.tps);
%computation of v11M
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
range = (lds.tps-1)*lds.nphase+(lds.phases);
jac(range,[lds.phases range]) = bordBVP_R2_bc1;

f1 = -jac(1:lds.ncoords-lds.nphase,lds.ncoords+1);
%computation of v1
b = []; b(lds.ncoords+1)=1; b=b';
jace=[jac(lds.coords,lds.coords),rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

jace=[jac(lds.coords,lds.coords) p;q' 0];
jace2 = jace;
vext = jace\b;
v1 = vext(lds.coords);
vps = reshape(v1,lds.nphase,lds.tps);
% rescale vext
ficd = dot(vps,vps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v1 = v1/sqrt(ic);

%computation of v1*
wext = jace'\b;
w1 = wext(1:lds.ncoords-lds.nphase)';

% computation of v2
v1ps = reshape(v1,lds.nphase,lds.tps);
ic = zeros(1,lds.ncoords);
range0 = lds.cols_p1;
range1 = lds.cols_p1_coords;
range2 = lds.cols;
range3 = lds.phases;
range6 = lds.col_coords;
wt = lds.wt';
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
for j=lds.tsts
    xp = ups(:,range0)*lds.wt;
    v1_1 = v1(range1);
    v1psC(:,range2) = v1ps(:,range0)*lds.wt;
    v1psW1(:,range2) = (v1psC(:,range2).*(repmat(gl_weight(lds.ncol)',lds.nphase,1)/2))*lds.dt(j); %compose v1W1   
    p = v1psW1(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(p,1,(lds.nphase*(lds.ncol+1)));    
    range5 = lds.phases; 
    for c=lds.cols
        xt = xp(:,c);
        stenv1v1 = zeros(lds.nphase,lds.nphase);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        sjac  = cjac(lds.func,lds.Jacobian,xt,pt,lds.ActiveParams);
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);        
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        sysjac(range5,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',sjac);
        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            for d2=lds.phases
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
            end            
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
        end
        fxhessv1v1(range3)=wtk.*shv1(:,kr2)*v1_1;%B_CxMxMv1Mv1M
        fxtensv1v1v1(range3)=wtk.*stenv1v1(:,kr2)*v1_1;
        range3 = range3+lds.nphase;
        range5 = range5+lds.nphase;
    end
    fxjac(range6,range1)=sysjac;
    
    range0 = range0 + lds.ncol;
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;    
    range6 = range6 + lds.ncol_coord;
end


%add integral constraint
jace(lds.ncoords+1,[lds.coords]) = ic;
v1C = reshape(v1psC,lds.ncoords-lds.nphase,1);

% fprintf('computation v2')
rl = [];
rl(1:lds.ncoords-lds.nphase) = -T*v1C; rl(lds.ncoords+1) = 0;
v2 = jace\rl';
v2 = v2(1:lds.ncoords);
v2ps = reshape(v2,lds.nphase,lds.tps);

range1 = lds.cols_p1;
range2 = lds.cols;
for j=lds.tsts
    v2psC(:,range2) = v2ps(:,range1)*lds.wt; %value in collocation points
    range1 = range1 + lds.ncol;
    range2 = range2 + lds.ncol;
end

%rescale v1*
v2C = reshape(v2psC,lds.ncoords-lds.nphase,1);
ic = w1*v2C;
w1 = w1'/ic;

w1ps = reshape(w1,lds.nphase,lds.tps-1);

% change boundary conditions
range = (lds.tps-1)*lds.nphase + (lds.phases);
jace = jac;
jace(range,[lds.phases range]) = BVP_LC1_jac_bc;
J=[jace(lds.coords,lds.coords),rand(lds.ncoords,1);rand(1,lds.ncoords),0];
b = []; b(lds.ncoords+1)=1; b=b';
q=J\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=J'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);
jace = [jace(lds.coords,lds.coords) p;q' 0];

%computation of phi*
phi = jace'\b;
phi = phi(1:lds.ncoords-lds.nphase);
tmp = phi'*f1;
phi = phi/tmp;

% fprintf('computation of alpha1')
alpha1 = 1/2*phi'*fxhessv1v1';

%computation of h20
%add integral constraint
phips = reshape(phi,lds.nphase,lds.tps-1);
ic1 = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
wt = lds.wt';
for j=lds.tsts
    p1 = phips(:,range2)*wt;
    ic1(range1) = ic1(range1)+reshape(p1,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end


range0 = lds.cols_p1;
range1 = lds.cols_p1_coords;
range3 = lds.phases;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
for j=lds.tsts
    xp = ups(:,range0)*lds.wt;
    v1_1 = v1(range1);
    v2_1 = v2(range1);
    for c=lds.cols
        xt = xp(:,c);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);        
        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shv2(:,d1) = (wtk.*hess(:,kr2,d1))*v2_1;
        end
        fxhessv1v2(range3)=wtk.*shv1(:,kr2)*v2_1;%B_CxMxMv1Mv2M
        fxhessv2v2(range3)=wtk.*shv2(:,kr2)*v2_1;%B_CxMxMv2Mv2M
        range3 = range3+lds.nphase;
    end    
    range0 = range0 + lds.ncol;
    range1 = range1 + lds.ncol_coord;
end

% add integral constraint h20
jace(lds.ncoords+1,[lds.coords]) = ic1;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v1'-2*alpha1*T*f1; rl(lds.ncoords+1) = phi'*fxhessv1v2';
h20 = jace\rl';
h20 = h20(lds.coords);
h20ps = reshape(h20,lds.nphase,lds.tps);

% function
range1 = lds.cols_p1;
range2 = lds.phases;
range3 = lds.cols;
range4 = lds.cols_p1_coords;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
for j=lds.tsts
    h20psC(:,range3) = h20ps(:,range1)*lds.wt;
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;
    v1_1 = v1(range4);
    h20_1 = h20(range4);
    v2_1 = v2(range4);
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);       
        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
        end
        fxhessv1h20(range2)=wtk.*shv1(:,kr2)*h20_1;
        range2 = range2+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range3 = range3+lds.ncol;
    range4 = range4 + lds.ncol_coord;
end
h20C = reshape(h20psC,lds.ncoords-lds.nphase,1);

a1 = 1/6*w1'*(fxtensv1v1v1'+3*fxhessv1h20'-6*alpha1*fxjac*v1);

%computation of h11
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v2'-T*h20C; rl(lds.ncoords+1) = 1/2*phi'*fxhessv2v2';
h11 = jace\rl';
h11 = h11(lds.coords);

range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
range5 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
w1ps = reshape(w1,lds.nphase,lds.tps-1);
wt = lds.wt';
ic = zeros(1,lds.ncoords);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
for j=lds.tsts
    p = w1ps(:,range5)*wt;
    ic(range4) = ic(range4) + reshape(p,1,(lds.nphase*(lds.ncol+1)));
    % value of polynomial on each collocation point    
    xp = ups(:,range1)*lds.wt;
    v1_1 = v1(range4);
    v2_1 = v2(range4);
    h20_1 = h20(range4);
    h11_1 = h11(range4);
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        stenv1v1 = zeros(lds.nphase,lds.nphase);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);

        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shh20(:,d1) = (wtk.*hess(:,kr2,d1))*h20_1;
            for d2=lds.phases
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
            end    
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
        end
        fxhessv1h11(range2)=wtk.*shv1(:,kr2)*h11_1;
        fxhessh20v2(range2)=wtk.*shh20(:,kr2)*v2_1;
        fxtensv1v1v2(range2)=wtk.*stenv1v1(:,kr2)*v2_1;
        range2 = range2+lds.nphase;
    end   
    range1 = range1 + lds.ncol;
    range4 = range4 + lds.ncol_coord;
    range5 = range5 + lds.ncol;
end
rl = -T*ic;
rl(lds.ncoords+1) = 0;
jace2(1:lds.ncoords-lds.nphase,lds.ncoords+1) = v2C;
jace2(lds.ncoords-lds.nphase+1:lds.ncoords+1,lds.ncoords+1) = zeros(lds.nphase+1,1);
w2 = jace2'\rl';
w2 = w2(1:lds.ncoords-lds.nphase);

b = 1/(2*T)*w1'*(-2*alpha1*fxjac*v2+fxtensv1v1v2'+fxhessh20v2'+2*fxhessv1h11')...
    +1/(2*T)*w2'*(fxtensv1v1v1'+3*fxhessv1h20'-6*alpha1*fxjac*v1);

result = [a1/T b];
% ---------------------------------------------------------------
function [x,p,T] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global lds

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+1:lds.PeriodIdx+2);
x = x0(lds.coords);
T = x0(lds.PeriodIdx);
