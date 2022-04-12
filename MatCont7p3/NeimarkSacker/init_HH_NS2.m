function [x0,v0] = init_HH_NS2(odefile, x, p, s, ap, ntst, ncol, eps)
%
% [x0,v0] = init_HH_NS2(odefile, x, p, s, ap, ntst, ncol,eps)
%
global lds cds
% check input
if(size(ap)~= 2)
  error('Two active parameters are needed for a torus bifurcation curve continuation');
end
cds.curve = @equilibrium;
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_jacobian = curvehandles{4};
cds.curve_hessians = curvehandles{5};

x0 = x(:,s.index);

% initialize lds
init_lds(odefile,x,p,s,ap,ntst,ncol);

func_handles = feval(odefile);
symord = 0; 
symordp = 0;

if     ~isempty(func_handles{9}),   symord = 5; 
elseif ~isempty(func_handles{8}),   symord = 4; 
elseif ~isempty(func_handles{7}),   symord = 3; 
elseif ~isempty(func_handles{5}),   symord = 2; 
elseif ~isempty(func_handles{3}),   symord = 1; 
end
if     ~isempty(func_handles{6}),   symordp = 2; 
elseif ~isempty(func_handles{4}),   symordp = 1; 
end
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.symjac = 1;
cds.symhess = 0;


lds.odefile = odefile;
lds.func = func_handles{2};
lds.Jacobian  = func_handles{3};
lds.JacobianP = func_handles{4};
lds.Hessians  = func_handles{5};
lds.HessiansP = func_handles{6};
lds.Der3 = func_handles{7};
lds.Der4 = func_handles{8};
lds.Der5 = func_handles{9};

x = x0(1:lds.nphase);
p(ap) = x0(lds.nphase+1:lds.nphase+2);

eds.ActiveParams = ap;
eds.P0 = p;
lds.P0 = p;
cds.oldJac = [];
cds.oldJacX = [];

%Preparing eigenvectors and derivatives.
  pp=p;
  p = num2cell(p);
  nphase=lds.nphase;
  jac = cjac(lds.func,lds.Jacobian,x,p,lds.ActiveParams);
  [X,D] = eig(jac);
  D=diag(D);
  index=find(abs(real(D))<1e-6 & sign(imag(D))==1);			% This should give a 1x2 vector
  if(size(index)~=2)                                  	% Otherwise a neutral saddle might be involved.
    debug('Neutral saddle ?\n');
    x0 =[];
    v0 =[];
    return;
  end
  if(imag(D(index(1))) < imag(D(index(2)))) 			% swap if necessary so that omega1>omega2
    index = [index(2);index(1)];
  end
  ev1 = D(index(1));
  ev2 = D(index(2));
  q1 = X(:,index(1));
  q2 = X(:,index(2));
  [XX,DD] = eig(jac');
  DD=diag(DD);
  index2=find(abs(real(DD))<1e-6 & sign(imag(DD))== -1);		% This should again give a 1x2 vector
  if(imag(DD(index2(1))) > imag(DD(index2(2))))	 		% swap if necessary
    index2 = [index2(2);index2(1)];
  end
  qad1 = XX(:,index2(1));
  qad2 = XX(:,index2(2));
  p1=qad1/(q1'*qad1);
  p2=qad2/(q2'*qad2);
  hessIncrement = (cds.options.Increment)^(3.0/4.0);
  ten3Increment = (cds.options.Increment)^(3.0/5.0);
  if (cds.options.SymDerivative >= 2)
    hess = chess(odefile,lds.Jacobian,lds.Hessians,x,p,lds.ActiveParams);
  else
    hess = [];
  end
  if (cds.options.SymDerivative >= 3)
    tens = ctens3(odefile,lds.Jacobian,lds.Hessians,lds.Der3,x,p,lds.ActiveParams);
  else
    tens = [];
  end
%Checking nondegeneracy
%2nd order vectors
  h2000 = (2*ev1*eye(nphase)-jac)\multilinear2(lds.func,hess,q1,q1,x,p,hessIncrement);           % (2iw_1-A)\B(q1,q1)
  h1100 = -real(jac\multilinear2(lds.func,hess,q1,conj(q1),x,p,hessIncrement));                  % -A\B(q1,conj(q1))
  h1010 = ((ev1+ev2)*eye(nphase)-jac)\multilinear2(lds.func,hess,q1,q2,x,p,hessIncrement);       % (i(w_1+w_2)-A)\B(q1,q2)
  h1001 = ((ev1-ev2)*eye(nphase)-jac)\multilinear2(lds.func,hess,q1,conj(q2),x,p,hessIncrement);	% (i(w_1-w_2)-A)\B(q1,conj(q2))
  h0020 = (2*ev2*eye(nphase)-jac)\multilinear2(lds.func,hess,q2,q2,x,p,hessIncrement);           % (2iw_2-A)\B(q2,q2)
  h0011 = -real(jac\multilinear2(lds.func,hess,q2,conj(q2),x,p,hessIncrement));                  % -A\B(q2,conj(q2))
%3rd order vectors
  h2100 = multilinear3(lds.func,tens,q1,q1,conj(q1),x,p,ten3Increment);              %  C(q1,q1,conj(q1))
  h2100 = h2100 + 2*multilinear2(lds.func,hess,q1,h1100,x,p,hessIncrement);          %+2B(h1100,q1)
  h2100 = h2100 + multilinear2(lds.func,hess,h2000,conj(q1),x,p,hessIncrement);      %+ B(h2000,conj(q1))    
  g2100 = real(p1'*h2100/2.0);%f2100 = imag(p1'*h2100/2.0);
  h1110 = multilinear3(lds.func,tens,q2,q1,conj(q1),x,p,ten3Increment);              %  C(q2,q1,conj(q1))
  h1110 = h1110 + multilinear2(lds.func,hess,h1100,q2,x,p,hessIncrement);            %+ B(q2,h1100)
  h1110 = h1110 + multilinear2(lds.func,hess,conj(h1001),q1,x,p,hessIncrement);      %+ B(q1,conj(h1001))
  h1110 = h1110 + multilinear2(lds.func,hess,h1010,conj(q1),x,p,hessIncrement);      %+ B(conj(q1),h1010)
  g1110 = real(p2'*h1110);%f1110 = imag(p2'*h1110);
  h1011 = multilinear3(lds.func,tens,q1,q2,conj(q2),x,p,ten3Increment);              %  C(q1,q2,conj(q2))
  h1011 = h1011 + multilinear2(lds.func,hess,h0011,q1,x,p,hessIncrement);            %+ B(q1,h0011)
  h1011 = h1011 + multilinear2(lds.func,hess,h1001,q2,x,p,hessIncrement);            %+ B(q2,h1001)
  h1011 = h1011 + multilinear2(lds.func,hess,h1010,conj(q2),x,p,hessIncrement);      %+ B(conj(q2),h1010)
  g1011 = real(p1'*h1011);f1011=imag(p1'*h1011);
  h0021 = multilinear3(lds.func,tens,q2,q2,conj(q2),x,p,ten3Increment);              %  C(q2,q2,conj(q2))
  h0021 = h0021 + 2*multilinear2(lds.func,hess,q2,h0011,x,p,hessIncrement);          %+2B(h0011,q2)
  h0021 = h0021 + multilinear2(lds.func,hess,h0020,conj(q2),x,p,hessIncrement);      %+ B(h0020,conj(q2))
  g0021 = real(p2'*h0021/2.0);f0021=imag(p2'*h0021/2.0);
  if(abs(g2100*g0021-g1011*g1110)<1e-12)
    error('Bifurcation is degenerate!');
  end
%Checking Transversality
  J1 = cjacp(lds.func,lds.JacobianP,x,p,ap);
  temp = -jac\J1;
  s1 = [1;0]; s2 = [0;1];
  hessp=chessp(lds.func,lds.Jacobian,lds.HessiansP,x,p,ap);
  test11 = hessp(:,:,1)*q1 + multilinear2(lds.func,hess,q1,temp*s1,x,p,hessIncrement);
  test12 = hessp(:,:,2)*q1 + multilinear2(lds.func,hess,q1,temp*s2,x,p,hessIncrement);
  test21 = hessp(:,:,1)*q2 + multilinear2(lds.func,hess,q2,temp*s1,x,p,hessIncrement);
  test22 = hessp(:,:,2)*q2 + multilinear2(lds.func,hess,q2,temp*s2,x,p,hessIncrement);
  AA = [ p1'*test11 p1'*test12; p2'*test21 p2'*test22 ];
  newbase = x + (h0011- temp*inv(real(AA))*[g1011 ; g0021])*eps^2;
  x0 = zeros(ntst*nphase+2,1);v0 = zeros(ntst*nphase+2,1);
  for i=1:lds.ncoords/lds.nphase
    zz = exp(sqrt(-1.0)*2*pi*i/((lds.ncoords/lds.nphase)-1));
    x0(((i-1)*nphase+1):((i-1)*nphase+nphase)) = newbase + real(2*q2*zz*eps + h0020*zz^2*eps^2);
    v0(((i-1)*nphase+1):((i-1)*nphase+nphase)) = real(2*q2*zz + 2*h0020*zz^2*eps)+ 2*(newbase-x)/eps;
  end
% generate a new mesh and interpolate
resc=norm(v0);
[x0,v0]=new_mesh(x0,v0,ntst,ncol);

% append period, parameters and k
domega1 = (-imag(AA(1,:))*inv(real(AA))*[g1011;g0021] + f1011);
domega2 = (-imag(AA(2,:))*inv(real(AA))*[g1011;g0021] + f0021);
omega1 = imag(ev1)+domega1*eps^2;
omega2 = imag(ev2)+domega2*eps^2;
k = cos(2*pi*omega1/omega2);
x0(lds.ncoords+1)= (2*pi/omega2);
x0(lds.ncoords+(2:3)) = pp(ap)-(real(AA)\[g1011;g0021])*eps^2;
x0(lds.ncoords+4) = k;
v0(lds.ncoords+1)= -4*pi*domega2*eps/imag(ev2)^2;
v0(lds.ncoords+(2:3)) = -2*(real(AA)\[g1011;g0021])*eps;
v0(lds.ncoords+4) = -4*pi*(domega1*imag(ev2)-domega2*imag(ev1))*eps/imag(ev2)^2*sin(2*pi*omega1/omega2);
v0(lds.ncoords+(1:4))=v0(lds.ncoords+(1:4))/resc;
v0=v0/norm(v0);

[x1,p2,T]=rearr(x0);p1 = num2cell(p2);

jac = spalloc(2*lds.ncoords-lds.nphase,2*lds.ncoords-lds.nphase,(2*lds.tps-1)*lds.nphase*3+(2*lds.ncol+1)*lds.nphase*2*lds.ntst);
ups = reshape(x1,lds.nphase,lds.tps);
% function

range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
for ns = 1:2
    range0 = lds.cols_p1;
    for j = lds.tsts
        xp = ups(:,range0)*lds.wt;
        jac(range1,range2) = bordBVP_NS_f(lds.func,xp,p1,T,j);
        range0 = range0 + lds.ncol;
        range1 = range1 + lds.ncol_coord;
        range2 = range2 + lds.ncol_coord;
    end
end
% boundary conditions
range  = 2*(lds.tps-1)*lds.nphase+ (lds.phases);
range2  = 2*(lds.ncoords-lds.nphase)+lds.phases;
range1 = lds.ncoords-lds.nphase+lds.phases;
jac(range,[lds.phases range1 range2]) = bordBVP_NS_bc1(k);


%compute borders

dimdoub=2*lds.ncoords-lds.nphase;
jace=[jac rand(dimdoub,2);rand(2,dimdoub) zeros(2,2)];
b=zeros(dimdoub+2,2);b(dimdoub+1,1)=1;b(dimdoub+2,2)=1;
q=jace\b;q=q(1:dimdoub,:);q=orth(q);
lds.NS_phi0 = q(:,1); 
lds.NS_phi1 = q(:,2);
p=jace'\b;p=p(1:dimdoub,:);p=orth(p);
lds.NS_psi0 = p(:,1); 
lds.NS_psi1 = p(:,2);

for i=1:lds.tps
    lds.upoldp(:,i) = T*feval(lds.func, 0, ups(:,i), p1{:});
end

%compute indices

A = NS_BVP_jac('BVP_LC_jac_f','BVP_LC_jac_bc','BVP_LC_jac_ic',x1,p2,T,3,2);

[Q,R] = qr(A');
Bord  = [jac lds.NS_psi0 lds.NS_psi1; [lds.NS_phi0';lds.NS_phi1'] zeros(2)];
bunit = zeros((2*lds.tps-1)*lds.nphase+2,2);
bunit(end-1:end,:) = eye(2);
sn = Bord\bunit;
st = Bord'\bunit;
v = sn(1:end-2,:)';
w = st(1:end-2,:)';
w1 = w(:,end-lds.nphase+1:end);
% calculate g'
ups = reshape(x1,lds.nphase,lds.tps);

range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
cv=[];
t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
gx = zeros(4,lds.ncoords+3);gk=[];
for ns = 1:2
      range3 = lds.cols_p1_coords;
      v1     = cv  ;
      range0 = lds.cols_p1;
      for tstpt = lds.tsts
          xp  = ups(:,range0)*lds.wt;
          cv = v(:,range2)';
          cw = w(:,range1);
          range = lds.phases;
          for c = lds.cols
              xt    = xp(:,c);
              sysj  = cjac(lds.func,lds.Jacobian,xt,p1,lds.ActiveParams);
              sysh  = chess(lds.func,lds.Jacobian,lds.Hessians,xt,p1,lds.ActiveParams);
              syshp = chessp(lds.func,lds.Jacobian,lds.HessiansP,xt,p1,lds.ActiveParams);
              wtk   = lds.wt(kr1,c(ones(1,lds.nphase)))';
              for d = lds.phases
                  sh1(:,d) = (wtk.*sysh(:,kr2,d))*cv(:,1);
                  sh2(:,d) = (wtk.*sysh(:,kr2,d))*cv(:,2);
              end      
              t11 = T* wtk.*sh1(:,kr2);
              t21 = T* wtk.*sh2(:,kr2);
              t12 = (wtk.*sysj(:,kr2))*cv(:,1);
              t22 = (wtk.*sysj(:,kr2))*cv(:,2);
              t13 = T* wtk.*syshp(:,kr2,1)* cv(:,1);
              t23 = T* wtk.*syshp(:,kr2,1)* cv(:,2);
              t14 = T* wtk.*syshp(:,kr2,2)* cv(:,1);
              t24 = T* wtk.*syshp(:,kr2,2)* cv(:,2);
              syshess1(range,:) = [t11 t12 t13 t14];      
              syshess2(range,:) = [t21 t22 t23 t24];      
              range = range + lds.nphase;    
          end
          gx(1,[range3 lds.ncoords+(1:3)]) = gx(1,[range3 lds.ncoords+(1:3)]) + cw(1,:)*syshess1;
          gx(2,[range3 lds.ncoords+(1:3)]) = gx(2,[range3 lds.ncoords+(1:3)]) + cw(1,:)*syshess2;
          gx(3,[range3 lds.ncoords+(1:3)]) = gx(3,[range3 lds.ncoords+(1:3)]) + cw(2,:)*syshess1;
          gx(4,[range3 lds.ncoords+(1:3)]) = gx(4,[range3 lds.ncoords+(1:3)]) + cw(2,:)*syshess2;
          range0 = range0 + lds.ncol;
          range1 = range1 + lds.ncol_coord;
          range2 = range2 + lds.ncol_coord;
          range3 = range3 + lds.ncol_coord;
      end
end  
gk(1,1) = 2*w1(1,:)*v1(end-lds.nphase+1:end,1);
gk(2,1) = 2*w1(1,:)*v1(end-lds.nphase+1:end,2);
gk(3,1) = 2*w1(2,:)*v1(end-lds.nphase+1:end,1);
gk(4,1) = 2*w1(2,:)*v1(end-lds.nphase+1:end,2);
B = [A ; gx gk]*Q;
Jres = B(2+lds.ncoords:end,2+lds.ncoords:end)';
[Q,R,E] = qr(full(Jres));
index = [1 1;1 2;2 1;2 2];
[I,J] = find(E(:,1:2));
lds.index1 = index(I(J(1)),:);
lds.index2 = index(I(J(2)),:);
%----------------------------------------------------------------------
function [x,p,T] = rearr(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x) and parameters (p)
global lds
nap = length(lds.ActiveParams);
lds.PeriodIdx = lds.ncoords+1;

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+(1:nap));

x = x0(lds.coords);
T = x0(lds.ncoords+1);

lds.T =T;

%-----------------------------------------------------------------
function init_lds(odefile,x,p,s,ap,ntst,ncol)
global lds
lds=[];
lds.odefile = odefile;
func_handles = feval(lds.odefile);
lds.func = func_handles{2};
lds.Jacobian  = func_handles{3};
lds.JacobianP = func_handles{4};
lds.Hessians  = func_handles{5};
lds.HessiansP = func_handles{6};
lds.Der3=func_handles{7};
lds.Der4=func_handles{8};
lds.Der5=func_handles{9};
siz = size(func_handles,2);
if siz > 9
    j=1;
    for k=10:siz
        lds.user{j}= func_handles{k};
        j=j+1;
    end
else lds.user=[];
end
lds.nphase = size(s.data.evec,2);
lds.ActiveParams = ap;
lds.P0 = p;
set_ntst_ncol(ntst,ncol,(0:ntst)/ntst);
lds.T = [];
lds.cols_p1 = 1:(lds.ncol+1);
lds.cols_p1_coords = 1:(lds.ncol+1)*lds.nphase;
lds.ncol_coord = lds.ncol*lds.nphase;
lds.col_coords = 1:lds.ncol*lds.nphase;
lds.pars = lds.ncoords+(1:4);
lds.phases = 1:lds.nphase;
lds.ntstcol = lds.ntst*lds.ncol;
lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwwt = kron(lds.wt',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);

lds.PD_psi = [];
lds.PD_phi = [];
lds.PD_new_phi = [];
lds.PD_new_psi = [];
lds.PD_switch = 0;

lds.BP_psi = [];
lds.BP_phi = [];
lds.BP_psi1 = [];
lds.BP_phi1 = [];
lds.BP_new_phi = [];
lds.BP_new_psi = [];
lds.BP_new_psi1 = [];
lds.BP_new_phi1 = [];
lds.BP_switch = 0;

lds.BPC_switch = 0;
lds.BPC_psi = [];
lds.BPC_phi1 = [];
lds.BPC_phi2 = [];

lds.LPC_phi = [];
lds.LPC_psi = [];
lds.LPC_new_phi = [];
lds.LPC_new_psi = [];
lds.LPC_switch = 0;

lds.NS_psi0 = [];
lds.NS_psi1 = [];
lds.NS_phi0 = [];
lds.NS_phi1 = [];
lds.NS1_new_phi = [];
lds.NS2_new_phi = [];
lds.NS1_new_psi = [];
lds.NS2_new_psi = [];
lds.NS_new_phi = [];
lds.NS_new_psi = [];
lds.NS_switch = 0;
lds.NS1_switch = 0;
lds.NS2_switch = 0;

lds.bialt_M1 = [];
lds.bialt_M2 = [];
lds.bialt_M3 = [];
lds.bialt_M4 = [];
lds.multipliers = nan;
lds.monodromy = [];
lds.multi_r1 = [];
lds.multi_r2 = [];
lds.BranchParam = lds.ActiveParams;
lds.ups = [];
lds.vps = [];
lds.BranchParams=[];
lds.tsts = 1:lds.ntst;
lds.cols = 1:lds.ncol;
