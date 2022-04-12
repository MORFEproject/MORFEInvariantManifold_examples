function result=nf_CPC(x)
%
% calculates CPC normal form coefficient
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
range0 = lds.cols_p1;         % 1:(ncol+1)
range1 = lds.col_coords;      % 1:ncol*nphase
range2 = lds.cols_p1_coords;  % 1:(ncol+1)*nphase
for j=lds.tsts   % 1:ntest
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
jac(range,lds.coords) = bordBVP_CPC_bc2(lds.func,ups,pt);

J = [jac(lds.coords,lds.coords) rand(lds.ncoords,1);rand(1,lds.ncoords) 0];
b=zeros(lds.ncoords+1,1);
b(lds.ncoords+1,1)=1;
p=J'\b;
q=J\b;
p=p(1:lds.ncoords,1);p=p/norm(p);
q=q(1:lds.ncoords,1);q=q/norm(q);
J = [jac(:,lds.coords) [p;0]];
rl = -T*jac(1:lds.ncoords-lds.nphase,lds.ncoords+1);
rl = [rl;zeros(lds.nphase,1);0];
v1 = J\rl;
v1 = v1(lds.coords);

f1 = -jac(1:lds.ncoords-lds.nphase,end);

jace = J;
jace(lds.ncoords+1,lds.coords) = q';
b=zeros(lds.ncoords+1,1);
b(lds.ncoords+1,1)=1;
wext = jace'\b; %wext*f = 0
w1 = wext(1:lds.ncoords-lds.nphase)';
v1ps = reshape(v1,lds.nphase,lds.tps);
range1 = lds.cols_p1;
range2 = lds.cols;
for j=lds.tsts
    % value of polynomial on each collocation point
    v1Cps(:,range2) = v1ps(:,range1)*lds.wt;
    range1 = range1 + lds.ncol;
    range2 = range2 + lds.ncol;
end
v1C = reshape(v1Cps,lds.ncoords-lds.nphase,1);
ic = w1*v1C;
w1 = w1/ic;
w1ps = reshape(w1,lds.nphase,lds.tps-1);

ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
wt = lds.wt';
for j=lds.tsts
    p = w1ps(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end

b = [T*ic 0];
jacee = jace;
jacee(:,end) = [v1C;zeros(lds.nphase+1,1)];
w2 = jacee'\b';
w2 = w2(1:lds.ncoords-lds.nphase);
% range2 = 1:lds.ncol;
% for j=lds.tsts
%     w1Cps(:,range2) = (w1ps(:,range2)./(repmat(gl_weight(lds.ncol)',lds.nphase,1)/2))/lds.dt(j);        
%     range2 = range2 + lds.ncol;
% end
% w1C = reshape(w1Cps,lds.ncoords-lds.nphase,1);

w2ps = reshape(w2,lds.nphase,lds.tps-1);
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
        fxhessv1v1(range2)=wtk.*shv1(:,kr2)*v1_1;
        fxtens3v1v1v1(range2)=wtk.*stenv1v1(:,kr2)*v1_1;
        range2 = range2+lds.nphase;
        range5 = range5+lds.nphase;
    end   
    fxjac(range3,range4)=sysjac;
    range1 = range1+lds.ncol;
    range3 = range3 + lds.ncol_coord;
    range4 = range4 + lds.ncol_coord;
end

% range2 = lds.cols;
% for j=lds.tsts
%     v1psW1(:,range2) = (v1Cps(:,range2).*(repmat(gl_weight(lds.ncol)',lds.nphase,1)/2))*lds.dt(j); %compose v1W1   
%     range2 = range2+lds.ncol;
% end
% v1W = reshape(v1psW1,lds.ncoords-lds.nphase,1);
% %alpha1 = 1/2*v1W'*(fxhessv1v1'+2*fxjac*v1)+1;
alpha1 = 0;

% Check that we are in a CPC point
b=1/2*w1*(fxhessv1v1'+2*fxjac*v1);
if abs(b)>1e-6
    fprintf('Warning: the quadratic coefficient of the normal form b=%d is not zero.',b);
    fprintf('\n');
end

%computation of h2
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
wt = lds.wt';
for j=lds.tsts
    p = w2ps(:,range2)*wt;
    ic(range1) = ic(range1)+reshape(p,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end

jace(end,lds.coords) = ic;
rl = [T*fxhessv1v1'+2*T*(fxjac*v1+f1)-2*alpha1*T*f1;zeros(lds.nphase+1,1)];
h2 = jace\rl;
h2 = h2(lds.coords);

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
    h2_1 = h2(range4);    
    % evaluate function value on each collocation point
    for c=lds.cols                    
        xt = xp(:,c);       
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);        
        for d1=lds.phases
            shh2(:,d1) = (wtk.*hess(:,kr2,d1))*h2_1;
        end
        fxhessh2v1(range2)=wtk.*shv1(:,kr2)*v1_1;
        range2 = range2+lds.nphase;
    end   
    range1 = range1+lds.ncol;
    range4 = range4 + lds.ncol_coord;
end

result = 1/6*w1*(-6*alpha1*fxjac*v1+3*fxjac*h2+3*fxhessv1v1'+6*fxjac*v1+3*fxhessh2v1'+fxtens3v1v1v1');

% ---------------------------------------------------------------
function [x,p,T] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global lds

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+1:lds.PeriodIdx+2);
x = x0(lds.coords);
T = x0(lds.PeriodIdx);
