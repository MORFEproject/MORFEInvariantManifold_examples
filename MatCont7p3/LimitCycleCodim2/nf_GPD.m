function result = nf_GPD(x)
%
% calculates generalized flip normal form coefficient
%

global lds
[xt,p,T] = rearr(x);
pt = num2cell(p);
ups = reshape(xt,lds.nphase,lds.tps);

J = BVP_jac('BVP_PD_jac',xt,p,T,1,1); 
[LJ,UJ] = lu(J);
b = []; b(lds.ncoords+1)=1; b=b';
wext = LJ'\(UJ'\b);
vext = UJ\(LJ\b);
% %compute borders
v1 = vext(1:lds.ncoords);
v1ps = reshape(v1,lds.nphase,lds.tps);
% rescale v1
ficd = dot(v1ps,v1ps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v1 = v1/sqrt(ic);
v1ps = reshape(v1,lds.nphase,lds.tps);

% rescale v1*
w1 = wext(1:lds.ncoords-lds.nphase);
range1 = lds.cols_p1;
range2 = lds.phases;
range3 = lds.cols;
for j = lds.tsts
    xp = ups(:,range1)*lds.wt;   
    v1psC(:,range3) = v1ps(:,range1)*lds.wt;
    for c = lds.cols
        xt = xp(:,c);
        fC(range2) = feval(lds.func,0,xt,pt{:});
        range2 = range2 + lds.nphase;
    end    
    range1 = range1 + lds.ncol; 
    range3 = range3 + lds.ncol;
end
fC = fC';
v1C = reshape(v1psC,lds.ncoords-lds.nphase,1);
ic = w1'*v1C;
w1 = w1/ic;

%computation of phi*

% change boundary conditions
range = (lds.tps-1)*lds.nphase + (lds.phases);
jac = J;
jac(range,[lds.phases range]) = BVP_LC1_jac_bc;
% remove borders
jac(:,end)=[];jac(end,:)=[];

%compute borders
jace=[jac,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

%add borders
jac=[jac p;q' 0];
phi = jac'\b;
phi = phi(1:lds.ncoords-lds.nphase);
ic = phi'*fC;
phi = phi/ic;

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
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
        end
        fxhessv1v1(range2)=wtk.*shv1(:,kr2)*v1_1;
        range2 = range2+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
end

%computation of alpha1
alpha1 = 1/2*phi'*fxhessv1v1';

%computation of h22
%computation of integral constraint for h22 and h33
phips = reshape(phi,lds.nphase,lds.tps-1);
w1ps = reshape(w1,lds.nphase,lds.tps-1);
ic2 = zeros(1,lds.ncoords);
ic3 = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
wt = lds.wt';
for j=lds.tsts
    p =phips(:,range2)*wt;
    ic2(range1) = ic2(range1)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    p = w1ps(:,range2)*wt;
    ic3(range1) = ic3(range1)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end
%add integral constraint h22
jac(lds.ncoords+1,[lds.coords]) =ic2;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v1-2*alpha1*T*fC'; rl(lds.ncoords+1) = 0;
h22 = jac\rl';
h22 = h22(1:lds.ncoords);


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
    h22_1 = h22(range4);
    range5 = lds.phases;         
    % evaluate function value on each collocation point
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
        fxhessv1h22(range2)=wtk.*shv1(:,kr2)*h22_1;
        fxtensv1v1v1(range2)=wtk.*stenv1v1(:,kr2)*v1_1;
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    fxjac(range3,range4)=sysjac;
    range1 = range1+lds.ncol;
    range3 = range3 + lds.ncol_coord;
    range4 = range4 + lds.ncol_coord;
end

% Check that we are in a GPD point
cc=1/(3*T)*w1'*(fxtensv1v1v1+3*fxhessv1h22-6*alpha1*(fxjac*v1)')';
%if abs(cc)>1e-6
%    fprintf('Warning: the cubic coefficient of the normal form c=%d is not zero.',cc);  
%    fprintf('\n');
%end

%computation of h33
%add integral constraint h33
J(lds.ncoords+1,[lds.coords]) = ic3;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v1v1+3*T*fxhessv1h22-6*alpha1*T*(fxjac*v1)'; 
rl(lds.ncoords+1) = 0;
h33 = J\rl';
h33 = h33(1:lds.ncoords);

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
    h22_1 = h22(range4);
    h33_1 = h33(range4);
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        stenv1v1 = zeros(lds.nphase,lds.nphase);
        stens4v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens4v1v1v1 = zeros(lds.nphase,lds.nphase);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        tens4 = ctens4(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,xt,pt,lds.ActiveParams);

        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shh22(:,d1) = (wtk.*hess(:,kr2,d1))*h22_1;
            for d2=lds.phases
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
                for d3=lds.phases
                    stens4v1(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v1_1;
                end
                stens4v1v1(:,d2,d1) = stens4v1v1(:,d2,d1)+(wtk.*stens4v1(:,kr2,d2,d1))*v1_1;
            end  
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
            stens4v1v1v1(:,d1) = stens4v1v1v1(:,d1)+(wtk.*stens4v1v1(:,kr2,d1))*v1_1;
        end
        fxhessv1h33(range2)=wtk.*shv1(:,kr2)*h33_1;
        fxhessh22h22(range2)=wtk.*shh22(:,kr2)*h22_1;
        fxtensv1v1h22(range2)=wtk.*stenv1v1(:,kr2)*h22_1;
        fxtens4v1v1v1v1(range2)=wtk.*stens4v1v1v1(:,kr2)*v1_1;
        range2 = range2+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
end
%computation of alpha2
alpha2 = 1/24*phi'*(fxtens4v1v1v1v1'+6*fxtensv1v1h22'+3*fxhessh22h22'+4*fxhessv1h33'...
    -12*alpha1*fxjac*h22-12*alpha1*fxhessv1v1')+alpha1^2;

%computation of h44
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v1v1+6*T*fxtensv1v1h22+3*T*fxhessh22h22...
    +4*T*fxhessv1h33-12*alpha1*T*((fxjac*h22)'+fxhessv1v1-2*alpha1*fC')-24*alpha2*T*fC';
rl(lds.ncoords+1) = 0;
h44 = jac\rl';
h44 = h44(lds.coords);

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
    h22_1 = h22(range4);
    h33_1 = h33(range4);
    h44_1 = h44(range4);
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);
        stens5v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase,lds.nphase);
        sten4v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stens5v1v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase); 
        stenv1v1 = zeros(lds.nphase,lds.nphase); 
        stenv1h22 = zeros(lds.nphase,lds.nphase);
        sten4v1v1v1 = zeros(lds.nphase,lds.nphase);
        stens5v1v1v1v1 = zeros(lds.nphase,lds.nphase);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        tens4 = ctens4(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,xt,pt,lds.ActiveParams);
        tens5 = ctens5(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,lds.Der5,xt,pt,lds.ActiveParams);

        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shh22(:,d1) = (wtk.*hess(:,kr2,d1))*h22_1;
            for d2=lds.phases                
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
                for d3=lds.phases
                    stens4v1(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v1_1;
                    for d4=lds.phases
                        stens5v1(:,d4,d3,d2,d1) = (wtk.*tens5(:,kr2,d1,d2,d3,d4))*v1_1;
                    end
                    stens5v1v1(:,d3,d2,d1) = stens5v1v1(:,d3,d2,d1)+wtk.*stens5v1(:,kr2,d3,d2,d1)*v1_1;
                end                
                sten4v1v1(:,d2,d1)=sten4v1v1(:,d2,d1)+wtk.*stens4v1(:,kr2,d2,d1)*v1_1;                
                stens5v1v1v1(:,d2,d1) = stens5v1v1v1(:,d2,d1)+wtk.*stens5v1v1(:,kr2,d2,d1)*v1_1;
            end  
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
            stenv1h22(:,d1)=stenv1h22(:,d1)+wtk.*stensv1(:,kr2,d1)*h22_1;
            sten4v1v1v1(:,d1)=sten4v1v1v1(:,d1)+wtk.*sten4v1v1(:,kr2,d1)*v1_1;                
            stens5v1v1v1v1(:,d1) = stens5v1v1v1v1(:,d1)+wtk.*stens5v1v1v1(:,kr2,d1)*v1_1;
        end
        fxhessv1h44(range2)=wtk.*shv1(:,kr2)*h44_1;
        fxhessh22h33(range2)=wtk.*shh22(:,kr2)*h33_1;
        fxtensv1v1h33(range2)=wtk.*stenv1v1(:,kr2)*h33_1;
        fxtensv1h22h22(range2)=wtk.*stenv1h22(:,kr2)*h22_1;
        fxtens4v1v1v1h22(range2)=wtk.*sten4v1v1v1(:,kr2)*h22_1;
        fxtens5v1v1v1v1v1(range2)=wtk.*stens5v1v1v1v1(:,kr2)*v1_1;        
        range2 = range2+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
end
%computation of e
result = 1/(120*T^2)*w1'*(fxtens5v1v1v1v1v1'+10*fxtens4v1v1v1h22'+15*fxtensv1h22h22'...
    +10*fxtensv1v1h33'+10*fxhessh22h33'+5*fxhessv1h44'-120*alpha2*fxjac*v1-20*alpha1*fxjac*h33);

%% ---------------------------------------------------------------
function [x,p,T] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global lds

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+1:lds.PeriodIdx+2);
x = x0(lds.coords);
T = x0(lds.PeriodIdx);

% -------------------------------------------------------------

function jac = BVP_jac(BVP_func,x,p,T,pars,nc)
global lds 

p2 = num2cell(p);
jac = feval(BVP_func,lds.func,x,p,T,pars,nc,lds,p2,lds.Jacobian,lds.ActiveParams,lds.JacobianP); 

