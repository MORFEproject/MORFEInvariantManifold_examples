function result = nf_CH(x0)
%
% calculates CH normal form coefficient
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

%computing vext and wext
[Jjac,Jjac2] = BVP_jac('bordBVP_R2_f','BVP_LC1_jac_bc',xt,p,T,theta,1,1);

b = []; b(lds.ncoords+1)=1; b=b';

jace=[Jjac,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);

Jjac=[Jjac p;q' 0];
[LJ,UJ] = lu(Jjac);

vext = UJ \(LJ\b);
ups = reshape(xt,lds.nphase,lds.tps);
%rescale vext
v1 = vext(1:lds.ncoords);
vps = reshape(v1,lds.nphase,lds.tps);
ficd = dot(vps,vps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v1 = v1/sqrt(ic);

wext = LJ'\(UJ'\b);
w1 = wext(1:lds.ncoords-lds.nphase);

% function
range1 = lds.cols_p1;
range2 = lds.phases;
range3 = lds.col_coords;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
v1ps = reshape(v1,lds.nphase,lds.tps);
for jk=lds.tsts
    % value of polynomial on each collocation point
    xp(:,range6) = ups(:,range1)*lds.wt;
    xpp = xp(:,range6);
    v1psC(:,range6) = v1ps(:,range1)*lds.wt;
    v1_1 = v1(range4);
    v11_1 = conj(v1_1);
    range5 = lds.phases;         
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xpp(:,c);
        stenv1v1 = zeros(lds.nphase,lds.nphase);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        sjac  = cjac(lds.func,lds.Jacobian,xt,pt,lds.ActiveParams);
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        f1(range2) = feval(lds.func, 0,  xt, pt{:});
        sysjac(range5,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',sjac);        
        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            for d2=lds.phases
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
            end    
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
        end
        fxhessv1v1(range2) = wtk.*shv1(:,kr2)*v1_1;
        fxhessv1v11(range2)= wtk.*shv1(:,kr2)*v11_1;
        fxtensv1v1v1(range2) = wtk.*stenv1v1(:,kr2)*v1_1;
        fxtensv1v1v11(range2) = wtk.*stenv1v1(:,kr2)*v11_1;
        range2 = range2 + lds.nphase;
        range5 = range5 + lds.nphase;
    end   
    fxjac(range3,range4) = sysjac;
    range1 = range1 + lds.ncol;
    range3 = range3 + lds.ncol_coord;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end
v1C = reshape(v1psC,lds.ncoords-lds.nphase,1);
f1 = f1';

ic = w1'*v1C;
w1 = w1/conj(ic);

%calculate phi*
jace=[Jjac2,rand(lds.ncoords,1);rand(1,lds.ncoords),0];
q=jace\b;
q=q(1:lds.ncoords,1);q=q/norm(q);
p=jace'\b;
p=p(1:lds.ncoords,1);p=p/norm(p);
Jjac2=[Jjac2 p;q' 0];
phi = Jjac2'\b;
%rescale phi
phi = phi(1:lds.ncoords-lds.nphase)';
phi = phi/(phi*f1);

%calculate h20
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
range3 = lds.cols;
twitheta = 2*1i*theta*lds.pwwt;
J = Jjac;
for jk  = lds.tsts
  % evaluate part of Jacobian matrix
  J(range1,range2) = bordBVP_R2_f(lds.func,xp(:,range3),pt,T,jk)+twitheta;  
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
  range3 = range3 + lds.ncol;
end

% remove borders
J(:,end)=[];J(end,:)=[];
rl=[];
rl(1:lds.ncoords-lds.nphase) = T*fxhessv1v1; rl(lds.ncoords) = 0;
h20 =J\rl.';

%computation of alpha1
alpha1 = phi*real(fxhessv1v11');

%computation of h11
phips = reshape(phi,lds.nphase,lds.tps-1);
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = lds.cols;
wt = lds.wt';
for jk=lds.tsts
    pw = phips(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(pw,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end

%add integral constraint h11
Jjac2(lds.ncoords+1,lds.coords) = ic;
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxhessv1v11')-alpha1*T*f1; rl(lds.ncoords+1) = 0;
h11 = Jjac2\rl';
h11 = h11(lds.coords);

range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
for jk=lds.tsts
    % value of polynomial on each collocation point
    xpp = ups(:,range1)*lds.wt;
    v1_1 = v1(range4);
    v11_1 = conj(v1_1);
    h11_1 = h11(range4);
    h20_1 = h20(range4);
    % evaluate function value on each collocation point
    for c=lds.cols            
        xt = xpp(:,c);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams); 
        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;            
            shv11(:,d1) = (wtk.*hess(:,kr2,d1))*v11_1;
            shh20(:,d1) = (wtk.*hess(:,kr2,d1))*h20_1;            
        end
        fxhessv1h11(range2) = wtk.*shv1(:,kr2)*h11_1;
        fxhessv1h20(range2) = wtk.*shv1(:,kr2)*h20_1;
        fxhessv11h20(range2) = wtk.*shv11(:,kr2)*h20_1;
        fxhessh20v11(range2) = wtk.*shh20(:,kr2)*v11_1;
        range2 = range2+lds.nphase;
    end   
    range1 = range1 + lds.ncol;
    range4 = range4 + lds.ncol_coord;
end

%computation of c1
c1 = -i/2*w1'*(fxtensv1v1v11.'+2*fxhessv1h11.'+fxhessv11h20.'-2*alpha1*fxjac*v1)+alpha1*theta/T;

% Check that we are in a CH point
if abs(real(c1))>1e-6
    fprintf('Warning: the real part of the cubic coefficient of the normal form %d is not close to zero.',real(c1));
    fprintf('\n');
end

%computation of h21
w1ps = reshape(w1,lds.nphase,lds.tps-1);
w11ps = conj(w1ps);
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = lds.cols;
wt = lds.wt';
for jk=lds.tsts
    pw = w11ps(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(pw,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end

%add integral constraint h21
Jjac(lds.ncoords+1,[lds.coords]) = ic;
rl = [];
rl(1:lds.ncoords-lds.nphase) = 2*T*fxhessv1h11.'+T*fxtensv1v1v11.'+T*fxhessh20v11.'-2*i*c1*T*v1C-2*alpha1*(T*fxjac*v1-i*theta*v1C);
rl(lds.ncoords+1) = 0;

h21 = Jjac\rl.';
h21 = h21(lds.coords);

%computation of h30
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
range3 = lds.cols;
twitheta = 3*1i*theta*lds.pwwt;
for jk  = lds.tsts
  % evaluate part of Jacobian matrix
  Jjac(range1,range2) = bordBVP_PD_jac_f(lds.func,xp(:,range3),pt,T,jk)+twitheta;  
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
  range3 = range3 + lds.ncol;
end

% remove borders
Jjac(:,end)=[];Jjac(end,:)=[];
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtensv1v1v1.'+3*T*fxhessv1h20.';
rl(lds.ncoords) = 0;
h30 = Jjac\rl.';

%computation of h31
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
    xpp = ups(:,range1)*lds.wt;
    h20Cps(:,range6) = h20ps(:,range1)*lds.wt;
    v1_1 = v1(range4);
    v11_1 = conj(v1_1);
    h11_1 = h11(range4);
    h20_1 = h20(range4);
    h02_1 = conj(h20_1);
    h21_1 = h21(range4);
    h12_1 = conj(h21_1);
    h30_1 = h30(range4);
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xpp(:,c);
        sten4v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase);
        stenv1v1 = zeros(lds.nphase,lds.nphase); 
        stenv1v11 = zeros(lds.nphase,lds.nphase);
        stenv11v11 = zeros(lds.nphase,lds.nphase);
        sten4v1v1v1 = zeros(lds.nphase,lds.nphase);
        sten4v1v1v11 = zeros(lds.nphase,lds.nphase);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        tens4 = ctens4(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,xt,pt,lds.ActiveParams);

        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shv11(:,d1) = (wtk.*hess(:,kr2,d1))*v11_1;
            shh11(:,d1) = (wtk.*hess(:,kr2,d1))*h11_1;
            shh02(:,d1) = (wtk.*hess(:,kr2,d1))*h02_1;
            for d2=lds.phases                
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
                stensv11(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v11_1;
                for d3=lds.phases
                    stens4v1(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v1_1;
                end                
                sten4v1v1(:,d2,d1)=sten4v1v1(:,d2,d1)+wtk.*stens4v1(:,kr2,d2,d1)*v1_1;                
            end  
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
            stenv1v11(:,d1)=stenv1v11(:,d1)+wtk.*stensv1(:,kr2,d1)*v11_1;
            stenv11v11(:,d1)=stenv11v11(:,d1)+wtk.*stensv11(:,kr2,d1)*v11_1;
            sten4v1v1v1(:,d1)=sten4v1v1v1(:,d1)+wtk.*sten4v1v1(:,kr2,d1)*v1_1;                
            sten4v1v1v11(:,d1)=sten4v1v1v11(:,d1)+wtk.*sten4v1v1(:,kr2,d1)*v11_1;                
        end
        fxhessv1h21(range2)=wtk.*shv1(:,kr2)*h21_1;
        fxhessv1h12(range2)=wtk.*shv1(:,kr2)*h12_1;        
        fxhessv11h30(range2)=wtk.*shv11(:,kr2)*h30_1;
        fxhessh11h20(range2)=wtk.*shh11(:,kr2)*h20_1;        
        fxhessh11h11(range2)=wtk.*shh11(:,kr2)*h11_1;                
        fxhessv11h21(range2)=wtk.*shv11(:,kr2)*h21_1;
        fxhessh02h20(range2)=wtk.*shh02(:,kr2)*h20_1;
        fxtensv1v1h11(range2)=wtk.*stenv1v1(:,kr2)*h11_1;
        fxtensv1v1h02(range2)=wtk.*stenv1v1(:,kr2)*h02_1;
        fxtensv1v11h20(range2)=wtk.*stenv1v11(:,kr2)*h20_1;
        fxtensv1v11h11(range2)=wtk.*stenv1v11(:,kr2)*h11_1;
        fxtensv11v11h20(range2)=wtk.*stenv11v11(:,kr2)*h20_1;
        fxtens4v1v1v1v11(range2)=wtk.*sten4v1v1v1(:,kr2)*v11_1;
        fxtens4v1v1v11v11(range2)=wtk.*sten4v1v1v11(:,kr2)*v11_1;
        range2 = range2+lds.nphase;
    end   
    range1 = range1 + lds.ncol;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end
h20C = reshape(h20Cps,lds.ncoords-lds.nphase,1);

range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
range3 = lds.cols;
twitheta = 2*1i*theta*lds.pwwt;
for jk  = lds.tsts
  % evaluate part of Jacobian matrix
  Jjac(range1,range2) = bordBVP_PD_jac_f(lds.func,xp(:,range3),pt,T,jk)+twitheta;  
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
  range3 = range3 + lds.ncol;
end
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*fxtens4v1v1v1v11.'+3*T*fxtensv1v1h11.'+3*T*fxtensv1v11h20.'...
    +3*T*fxhessh11h20.'+3*T*fxhessv1h21.'+T*fxhessv11h30.'-6*i*c1*T*h20C-3*alpha1*(T*fxjac*h20-2*i*theta*h20C+T*fxhessv1v1.');
rl(lds.ncoords) = 0;
h31 = Jjac\rl.';
alpha2 = 1/4*phi*(real(fxtens4v1v1v11v11)'+fxtensv1v1h02.'+4*fxtensv1v11h11.'+2*fxhessh11h11.'...
    +2*fxhessv1h12.'+fxtensv11v11h20.'+fxhessh02h20.'+2*fxhessv11h21.'-4*alpha1*(fxjac*h11+real(fxhessv1v11)'))+alpha1^2;

%computation of h22
rl = [];
rl(1:lds.ncoords-lds.nphase) = T*real(fxtens4v1v1v11v11)'+T*fxtensv1v1h02.'-4*alpha2*T*f1...
    +4*T*real(fxtensv1v11h11)'+2*T*fxhessh11h11'+2*T*fxhessv1h12.'+T*fxtensv11v11h20.'...
    +T*real(fxhessh02h20)'+2*T*fxhessv11h21.'-4*alpha1*T*(fxjac*h11+real(fxhessv1v11)'-alpha1*f1);
rl(lds.ncoords+1) = 0;
h22 = Jjac2\rl.';
h22 = h22(lds.coords);

range1 = lds.cols_p1;
range2 = lds.phases;
range4 = lds.cols_p1_coords;
range6 = lds.cols;
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
h21ps = reshape(h21,lds.nphase,lds.tps);
for j=lds.tsts
    % value of polynomial on each collocation point
    xpp = ups(:,range1)*lds.wt;
    h21Cps(:,range6) = h21ps(:,range1)*lds.wt;
    v1_1 = v1(range4);
    v11_1 = conj(v1_1);
    h11_1 = h11(range4);
    h21_1 = h21(range4);
    h12_1 = conj(h21_1);
    h31_1 = h31(range4);
    h30_1 = h30(range4);
    h20_1 = h20(range4);
    h02_1 = conj(h20_1);
    h22_1 = h22(range4);
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xpp(:,c);

        stens5v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase,lds.nphase);
        sten4v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase);    
        sten4v1v11 = zeros(lds.nphase,lds.nphase,lds.nphase);        
        stens5v1v1v1 = zeros(lds.nphase,lds.nphase,lds.nphase);        
        stenv1v1 = zeros(lds.nphase,lds.nphase);         
        stenv1v11 = zeros(lds.nphase,lds.nphase);         
        stenv11v11 = zeros(lds.nphase,lds.nphase);         
        stenv1h02 = zeros(lds.nphase,lds.nphase);         
        stenv1h11 = zeros(lds.nphase,lds.nphase);         
        stenv11h11 = zeros(lds.nphase,lds.nphase);         
        sten4v1v1v1 = zeros(lds.nphase,lds.nphase);
        sten4v1v1v11 = zeros(lds.nphase,lds.nphase);        
        sten4v1v11v11 = zeros(lds.nphase,lds.nphase);        
        stens5v1v1v1v11 = zeros(lds.nphase,lds.nphase);

        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
        tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
        tens4 = ctens4(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,xt,pt,lds.ActiveParams);
        tens5 = ctens5(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,lds.Der4,lds.Der5,xt,pt,lds.ActiveParams);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        for d1=lds.phases
            shv1(:,d1) = (wtk.*hess(:,kr2,d1))*v1_1;
            shv11(:,d1) = (wtk.*hess(:,kr2,d1))*v11_1;
            shh11(:,d1) = (wtk.*hess(:,kr2,d1))*h11_1;
            shh12(:,d1) = (wtk.*hess(:,kr2,d1))*h12_1;
            shh20(:,d1) = (wtk.*hess(:,kr2,d1))*h20_1;
            shh02(:,d1) = (wtk.*hess(:,kr2,d1))*h02_1;
            for d2=lds.phases                
                stensv1(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v1_1;
                stensv11(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v11_1;
                for d3=lds.phases
                    stens4v1(:,d3,d2,d1) = (wtk.*tens4(:,kr2,d1,d2,d3))*v1_1;
                    for d4=lds.phases
                        stens5v1(:,d4,d3,d2,d1) = (wtk.*tens5(:,kr2,d1,d2,d3,d4))*v1_1;
                    end
                    stens5v1v1(:,d3,d2,d1) = stens5v1v1(:,d3,d2,d1)+wtk.*stens5v1(:,kr2,d3,d2,d1)*v1_1;
                end                
                sten4v1v1(:,d2,d1)=sten4v1v1(:,d2,d1)+wtk.*stens4v1(:,kr2,d2,d1)*v1_1;                
                sten4v1v11(:,d2,d1)=sten4v1v11(:,d2,d1)+wtk.*stens4v1(:,kr2,d2,d1)*v11_1;                
                stens5v1v1v1(:,d2,d1) = stens5v1v1v1(:,d2,d1)+wtk.*stens5v1v1(:,kr2,d2,d1)*v1_1;
            end  
            stenv1v1(:,d1)=stenv1v1(:,d1)+wtk.*stensv1(:,kr2,d1)*v1_1;
            stenv1v11(:,d1)=stenv1v11(:,d1)+wtk.*stensv1(:,kr2,d1)*v11_1;
            stenv11v11(:,d1)=stenv11v11(:,d1)+wtk.*stensv11(:,kr2,d1)*v11_1;
            stenv1h02(:,d1)=stenv1h02(:,d1)+wtk.*stensv1(:,kr2,d1)*h02_1;
            stenv1h11(:,d1)=stenv1h11(:,d1)+wtk.*stensv1(:,kr2,d1)*h11_1;
            stenv11h11(:,d1)=stenv11h11(:,d1)+wtk.*stensv11(:,kr2,d1)*h11_1;
            sten4v1v1v1(:,d1)=sten4v1v1v1(:,d1)+wtk.*sten4v1v1(:,kr2,d1)*v1_1;                
            sten4v1v1v11(:,d1)=sten4v1v1v11(:,d1)+wtk.*sten4v1v1(:,kr2,d1)*v11_1;                
            sten4v1v11v11(:,d1)=sten4v1v11v11(:,d1)+wtk.*sten4v1v11(:,kr2,d1)*v11_1;                
            stens5v1v1v1v11(:,d1) = stens5v1v1v1v11(:,d1)+wtk.*stens5v1v1v1(:,kr2,d1)*v11_1;
        end
        fxhessv1h11(range2)=wtk.*shv1(:,kr2)*h11_1;
        fxhessv1h22(range2)=wtk.*shv1(:,kr2)*h22_1;
        fxhessh11h21(range2)=wtk.*shh11(:,kr2)*h21_1;
        fxhessh12h20(range2)=wtk.*shh12(:,kr2)*h20_1;
        fxhessv11h31(range2)=wtk.*shv11(:,kr2)*h31_1;
        fxhessh20v11(range2)=wtk.*shh20(:,kr2)*v11_1;
        fxhessh02h30(range2)=wtk.*shh02(:,kr2)*h30_1;
        fxtensv1v1v11(range2)=wtk.*stenv1v1(:,kr2)*v11_1;
        fxtensv1v1h12(range2)=wtk.*stenv1v1(:,kr2)*h12_1;
        fxtensv1v11h21(range2)=wtk.*stenv1v11(:,kr2)*h21_1;
        fxtensv11v11h30(range2)=wtk.*stenv11v11(:,kr2)*h30_1;
        fxtensv1h02h20(range2)=wtk.*stenv1h02(:,kr2)*h20_1;
        fxtensv11h11h20(range2)=wtk.*stenv11h11(:,kr2)*h20_1;
        fxtensv1h11h11(range2)=wtk.*stenv1h11(:,kr2)*h11_1;
        fxtens4v1v1v1h02(range2)=wtk.*sten4v1v1v1(:,kr2)*h02_1;
        fxtens4v1v1v11h11(range2)=wtk.*sten4v1v1v11(:,kr2)*h11_1;
        fxtens4v1v11v11h20(range2)=wtk.*sten4v1v11v11(:,kr2)*h20_1;
        fxtens5v1v1v1v11v11(range2)=wtk.*stens5v1v1v1v11(:,kr2)*v11_1;        
        range2 = range2+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
    range6 = range6 + lds.ncol;
end
h21C = reshape(h21Cps,lds.ncoords-lds.nphase,1);

e1 = 1/(12*T^2)*w1'*(fxtens5v1v1v1v11v11.'+fxtens4v1v1v1h02.'+6*fxtens4v1v1v11h11.'...
    +6*fxtensv1h11h11.'+3*fxtensv1v1h12.'+3*fxtens4v1v11v11h20.'+3*fxtensv1h02h20.'...
    +6*fxtensv11h11h20.'+3*fxhessh12h20.'+6*fxtensv1v11h21.'+6*fxhessh11h21.'+3*fxhessv1h22.'...
    +fxtensv11v11h30.'+fxhessh02h30.'+2*fxhessv11h31.'-6*i*c1*h21C-12*alpha2*fxjac*v1...
    -6*alpha1*(fxjac*h21-i*theta/T*h21C+2*fxhessv1h11.'+fxtensv1v1v11.'+fxhessh20v11.'...
    -2*alpha1*fxjac*v1))+alpha2*i*theta/T^3+alpha1*i*c1/T^2-alpha1^2*i*theta/T^3;

result = real(e1);

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