function [x0,v0] = init_GH_LPC(odefile, x, p, s, ap, ntst, ncol, eps, varargin)
%
% [x0,v0] = init_GH_LPC(odefile, x, p, s, ap, ntst, ncol, eps)
%
global lds hds cds
% check input
if(size(ap)~= 2)
  error('Two active parameters are needed for a fold bifurcation curve continuation');
end
cds.curve = @equilibrium;
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_jacobian = curvehandles{4};
cds.curve_hessians = curvehandles{5};

x0 = x(:,s.index);

% initialize lds
if size(varargin,1)>0
    bp = varargin{1};
    init_lds(odefile,x,p,s,ap,ntst,ncol,bp);
else 
    init_lds(odefile,x,p,s,ap,ntst,ncol);
end

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

%Checking Nondegeneracy and the theoretical possibility to switch.
  pp=p;p = num2cell(p);nphase=size(x,1);
  jac = cjac(lds.func,lds.Jacobian,x,p,ap);
  [X,D] = eig(jac);
  index = find(abs(real(diag(D)))<1e-6 & sign(imag(diag(D)))==1);
if(isempty(index)==1)
  debug('Neutral saddle\n');
  return;
end
  ev0 = diag(D(index,index));
  q0 = X(:,index);
  q0 = q0/norm(q0);
  q1 = conj(q0);
  [X,D] = eig(jac');
  index = find(abs(real(diag(D)))<1e-6 & sign(imag(diag(D)))== -1);
  qad0 = X(:,index);
  p0 = qad0/(q0'*qad0);
  hessIncrement = (cds.options.Increment)^(3.0/4.0);
  ten3Increment = (cds.options.Increment)^(3.0/5.0);
  ten4Increment = (cds.options.Increment)^(3.0/6.0);
  ten5Increment = (cds.options.Increment)^(3.0/7.0);

  if (cds.options.SymDerivative >= 5)
      hess = chess(odefile,lds.Jacobian,lds.Hessians,x,p,lds.ActiveParams);
      tens = ctens3(odefile,lds.Jacobian,lds.Hessians,lds.Der3,x,p,lds.ActiveParams);      
      ten4 = feval(lds.Der4,0,x, p{:});
      ten5 = feval(lds.Der5,0,x, p{:});
  elseif (cds.options.SymDerivative >= 4)
      hess = chess(odefile,lds.Jacobian,lds.Hessians,x,p,lds.ActiveParams);
      tens = ctens3(odefile,lds.Jacobian,lds.Hessians,lds.Der3,x,p,lds.ActiveParams);      
      ten4 = feval(lds.Der4,0,x, p{:});
      ten5 = [];
  elseif (cds.options.SymDerivative >= 3)
      hess = chess(odefile,lds.Jacobian,lds.Hessians,x,p,lds.ActiveParams);
      tens = ctens3(odefile,lds.Jacobian,lds.Hessians,lds.Der3,x,p,lds.ActiveParams);      
      ten4 = [];
      ten5 = [];  
  elseif (cds.options.SymDerivative >= 2)
      hess = chess(odefile,lds.Jacobian,lds.Hessians,x,p,lds.ActiveParams);
      tens = [];      
      ten4 = [];
      ten5 = [];
  else
      hess = [];
      tens = [];      
      ten4 = [];
      ten5 = [];
  end

  Abor = [jac-ev0*eye(nphase) q0; p0' 0 ];
%-------------------------------------------------------------------------------
%2nd order vectors
  h20 = (2*ev0*eye(nphase)-jac)\multilinear2(lds.func,hess,q0,q0,x,p,hessIncrement); 
% (2iw-A)\B(q0,q0)
  h11 = -real(jac\multilinear2(lds.func,hess,q0,q1,x,p,hessIncrement)); 
% -A\B(q0,q1)
%3rd order vectors
  h30 = multilinear3(lds.func,tens,q0,q0,q0,x,p,ten3Increment); 
%  C(q0,q0,q1)
  h30 = h30 + 3*multilinear2(lds.func,hess,q0,h20,x,p,hessIncrement); 
%+3B(h20,q0)  
  h30 = (3*ev0*eye(nphase)-jac)\h30;
  h21 = multilinear3(lds.func,tens,q0,q0,q1,x,p,ten3Increment); 
%  C(q0,q0,q1)
  h21 = h21 + 2*multilinear2(lds.func,hess,q0,h11,x,p,hessIncrement); 
%+2B(h11,q0)
  h21 = h21 + multilinear2(lds.func,hess,h20,q1,x,p,hessIncrement); 
%+ B(h20,q1)
  g21 = p0'*h21/2.0;
  h21 = Abor\[ 2*g21*q0-h21 ; 0];
  h21 = h21(1:nphase);
%4th order vectors
  h31 = multilinear4(lds.func,ten4,q0,q0,q0,q1,x,p,ten4Increment); 
%  D(q0,q0,q0,q1)
  h31 = h31 + 3*multilinear3(lds.func,tens,q0,q0,h11,x,p,ten3Increment); 
%+3C(q0,q0,h11)
  h31 = h31 + 3*multilinear3(lds.func,tens,q0,q1,h20,x,p,ten3Increment); 
%+3C(q0,q1,h20)
  h31 = h31 + 3*multilinear2(lds.func,hess,h20,h11,x,p,hessIncrement); 
%+3B(h20,h11)
  h31 = h31 + 3*multilinear2(lds.func,hess,h21,q0,x,p,hessIncrement); 
%+3B(h21,q0)
  h31 = h31 +   multilinear2(lds.func,hess,h30,q1,x,p,hessIncrement); 
%+ B(h30,q1)
  h31 = (2*ev0*eye(nphase)-jac)\(h31 - 6*g21*h20);  
  h22 = multilinear4(lds.func,ten4,q0,q0,q1,q1,x,p,ten4Increment); 
%  D(q0,q0,q1,q1)
  h22 = h22 + 4*multilinear3(lds.func,tens,q0,q1,h11,x,p,ten3Increment); 
%+4C(q0,q1,h11)
  h22 = h22 + 2*real(multilinear3(lds.func,tens,q1,q1,h20,x,p,ten3Increment)); 
%+2*Re(C(q0,q0,h02))
  h22 = h22 + 4*real(multilinear2(lds.func,hess,h21,q1,x,p,hessIncrement)); 
%+2*Re(B(h21,q1))
  h22 = h22 + 2*multilinear2(lds.func,hess,h11,h11,x,p,hessIncrement); 
%+2B(h11,h11)
  h22 = -jac\(h22 + multilinear2(lds.func,hess,h20,conj(h20),x,p,hessIncrement)); %+ B(h20,h02)
%5th order rhs
  h32 = multilinear5(lds.func,ten5,q0,q0,q0,q1,q1,x,p,ten5Increment); 
%  E(q0,q0,q0,q1,q1)
  h32 = h32 + 6*multilinear4(lds.func,ten4,q0,q0,q1,h11,x,p,ten4Increment); 
%+6D(q0,q0,q1,h11)
  h32 = h32 + 3*multilinear4(lds.func,ten4,q1,q1,q0,h20,x,p,ten4Increment); 
%+3D(q1,q1,q0,h20)
  h32 = h32 + multilinear4(lds.func,ten4,q0,q0,q0,conj(h20),x,p,ten4Increment); %+ D(q0,q0,q0,h02)
  h32 = h32 + 6*multilinear3(lds.func,tens,h11,h11,q0,x,p,ten3Increment); 
%+6C(h11,h11,q0)
  h32 = h32 + 6*multilinear3(lds.func,tens,q1,h20,h11,x,p,ten3Increment); 
%+6C(q1,h20,h11)
  h32 = h32 + 6*multilinear3(lds.func,tens,q1,q0,h21,x,p,ten3Increment); 
%+6C(q1,q0,h21)
  h32 = h32 + 3*multilinear3(lds.func,tens,q0,h20,conj(h20),x,p,ten3Increment); 
%+3C(q1,h20,h02)
  h32 = h32 + 3*multilinear3(lds.func,tens,q0,q0,conj(h21),x,p,ten3Increment); 
%+3C(q0,q0,h12)
  h32 = h32 +  multilinear3(lds.func,tens,q1,q1,h30,x,p,ten3Increment); 
%+ C(q1,q1,h30)
  h32 = h32 + 6*multilinear2(lds.func,hess,h21,h11,x,p,hessIncrement); 
%+6B(h21,h11)
  h32 = h32 + 3*multilinear2(lds.func,hess,h22,q0,x,p,hessIncrement); 
%+3B(h22,q0)
  h32 = h32 + 3*multilinear2(lds.func,hess,h20,conj(h21),x,p,hessIncrement); 
%+3B(h12,h20)
  h32 = h32 + 2*multilinear2(lds.func,hess,h31,q1,x,p,hessIncrement); 
%+2B(h31,q1)
  h32 = h32 +   multilinear2(lds.func,hess,h30,conj(h20),x,p,hessIncrement); %+ B(h30,h02)
  lyap2 = real(p0'*h32/12);
%-------------------------------------------------------------------------------
%Checking Transversality : The new frame
  J1 = cjacp(lds.func,lds.JacobianP,x,p,ap);
  hessp=chessp(lds.func,lds.Jacobian,lds.HessiansP,x,p,ap);
  s1 = [1;0]; s2 = [0;1];
  h0010 = -jac\(J1*s1);
  h0001 = -jac\(J1*s2);
  test1 = hessp(:,:,1)*q0;
  test2 = hessp(:,:,2)*q0;
  test1 = test1 + multilinear2(lds.func,hess,q0,h0010,x,p,hessIncrement);
  test2 = test2 + multilinear2(lds.func,hess,q0,h0001,x,p,hessIncrement);
  gg1 = p0'*[test1 test2];
  hh1010 = Abor\[gg1(1)*q0-test1; 0];
  hh1001 = Abor\[gg1(2)*q0-test2; 0];
  h1010 = hh1010(1:nphase); h1001 = hh1001(1:nphase);
%-------------------------------------------------------------------------------
%Some derivatives wrt to parameters before we proceed:
%temp1 = B_1(q0,q0,s1)
%temp2 = B_1(q0,q1,s1)
%temp3 = C_1(q0,q0,q1,s1)+2B_1(h11,q0,s1)+B_1(h20,q1,s1)
%temp4 = B_1(q0,q0,s2)
%temp5 = B_1(q0,q1,s2)
%temp6 = C_1(q0,q0,q1,s2)+2B_1(h11,q0,s2)+B_1(h20,q1,s2)
%wrt to s1
  p1 = pp; p1(ap) = p1(ap) + cds.options.Increment*s1; p1=num2cell(p1);
  if (cds.options.SymDerivative >= 3)
    hess = chess(lds.func,lds.Jacobian,lds.Hessians,x,p1,lds.ActiveParams);
    tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,x,p1,lds.ActiveParams);
  end
  temp1 = multilinear2(lds.func,hess,q0,q0,x,p1,hessIncrement);
  temp2 = multilinear2(lds.func,hess,q0,q1,x,p1,hessIncrement);
  temp3 = multilinear3(lds.func,tens,q0,q0,q1,x,p1,ten3Increment);
  temp3 = temp3 + 2*multilinear2(lds.func,hess,h11,q0,x,p1,hessIncrement);
  temp3 = temp3 +   multilinear2(lds.func,hess,h20,q1,x,p1,hessIncrement);
  p1 = pp; p1(ap) = p1(ap) - cds.options.Increment*s1; p1=num2cell(p1);
  if (cds.options.SymDerivative >= 3)
    hess = chess(lds.func,lds.Jacobian,lds.Hessians,x,p1,lds.ActiveParams);
    tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,x,p1,lds.ActiveParams);
  end
  temp1 = temp1 - multilinear2(lds.func,hess,q0,q0,x,p1,hessIncrement);
  temp2 = temp2 - multilinear2(lds.func,hess,q0,q1,x,p1,hessIncrement);
  temp3 = temp3 - multilinear3(lds.func,tens,q0,q0,q1,x,p1,ten3Increment);
  temp3 = temp3 - 2*multilinear2(lds.func,hess,h11,q0,x,p1,hessIncrement);
  temp3 = temp3 -   multilinear2(lds.func,hess,h20,q1,x,p1,hessIncrement);
  temp1 = temp1/(2.0*cds.options.Increment);
  temp2 = temp2/(2.0*cds.options.Increment);
  temp3 = temp3/(2.0*cds.options.Increment);
%wrt to s2
  p1 = pp; p1(ap) = p1(ap) + cds.options.Increment*s2; p1=num2cell(p1);
  if (cds.options.SymDerivative >= 3)
    hess = chess(lds.func,lds.Jacobian,lds.Hessians,x,p1,lds.ActiveParams);
    tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,x,p1,lds.ActiveParams);
  end
  temp4 = multilinear2(lds.func,hess,q0,q0,x,p1,hessIncrement);
  temp5 = multilinear2(lds.func,hess,q0,q1,x,p1,hessIncrement);
  temp6 = multilinear3(lds.func,tens,q0,q0,q1,x,p1,ten3Increment);
  temp6 = temp6 + 2*multilinear2(lds.func,hess,h11,q0,x,p1,hessIncrement);
  temp6 = temp6 +   multilinear2(lds.func,hess,h20,q1,x,p1,hessIncrement);
  p1 = pp; p1(ap) = p1(ap) - cds.options.Increment*s2; p1=num2cell(p1);
  if (cds.options.SymDerivative >= 3)
    hess = chess(lds.func,lds.Jacobian,lds.Hessians,x,p1,lds.ActiveParams);
    tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,x,p1,lds.ActiveParams);
  end
  temp4 = temp4 - multilinear2(lds.func,hess,q0,q0,x,p1,hessIncrement);
  temp5 = temp5 - multilinear2(lds.func,hess,q0,q1,x,p1,hessIncrement);
  temp6 = temp6 - multilinear3(lds.func,tens,q0,q0,q1,x,p1,ten3Increment);
  temp6 = temp6 - 2*multilinear2(lds.func,hess,h11,q0,x,p1,hessIncrement);
  temp6 = temp6 -   multilinear2(lds.func,hess,h20,q1,x,p1,hessIncrement);
  temp4 = temp4/(2.0*cds.options.Increment);
  temp5 = temp5/(2.0*cds.options.Increment);
  temp6 = temp6/(2.0*cds.options.Increment);
  if (cds.options.SymDerivative >= 3)
    hess = chess(lds.func,lds.Jacobian,lds.Hessians,x,p,lds.ActiveParams);
    tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,x,p,lds.ActiveParams);
  end
%-------------------------------------------------------------------------------
%Expanding the Lyapunov Coefficient
  h2010 = temp1 + hessp(:,:,1)*h20 ;
  h2010 = h2010 + multilinear3(lds.func,tens,q0,q0,h0010,x,p,ten3Increment);
  h2010 = h2010 + multilinear2(lds.func,hess,h20,h0010,x,p,hessIncrement);
  h2010 = h2010 + 2*multilinear2(lds.func,hess,q0,h1010,x,p,hessIncrement);
  h2010 = (2*ev0*eye(nphase)-jac)\(h2010-2*h20*gg1(1));
  h1110 = temp2 + hessp(:,:,1)*h11 ;
  h1110 = h1110 + multilinear3(lds.func,tens,q0,q1,h0010,x,p,ten3Increment);
  h1110 = h1110 + multilinear2(lds.func,hess,h11,h0010,x,p,hessIncrement);
  h1110 = h1110 + 2*real(multilinear2(lds.func,hess,q1,h1010,x,p,hessIncrement));
  h1110 = -jac\(h1110-2*real(gg1(1))*h11);  
  h2110 = temp3 + hessp(:,:,1)*h21 ;
  h2110 = h2110 + multilinear4(lds.func,ten4,q0,q0,q1,h0010,x,p,ten4Increment);
  h2110 = h2110 + 2*multilinear3(lds.func,tens,q0,h11,h0010,x,p,ten3Increment);
  h2110 = h2110 + 2*multilinear3(lds.func,tens,q0,q1,h1010,x,p,ten3Increment);
  h2110 = h2110 + multilinear3(lds.func,tens,h20,q1,h0010,x,p,ten3Increment);
  h2110 = h2110 + multilinear3(lds.func,tens,q0,q0,conj(h1010),x,p,ten3Increment);
  h2110 = h2110 + 2*multilinear2(lds.func,hess,q0,h1110,x,p,hessIncrement);
  h2110 = h2110 + 2*multilinear2(lds.func,hess,h11,h1010,x,p,hessIncrement);
  h2110 = h2110 + multilinear2(lds.func,hess,h21,h0010,x,p,hessIncrement);
  h2110 = h2110 + multilinear2(lds.func,hess,h20,conj(h1010),x,p,hessIncrement);
  h2110 = h2110 + multilinear2(lds.func,hess,q1,h2010,x,p,hessIncrement);
% %
  h2001 = temp4 + hessp(:,:,2)*h20 ;
  h2001 = h2001 + multilinear3(lds.func,tens,q0,q0,h0001,x,p,ten3Increment);
  h2001 = h2001 + multilinear2(lds.func,hess,h20,h0001,x,p,hessIncrement);
  h2001 = h2001 + 2*multilinear2(lds.func,hess,q0,h1001,x,p,hessIncrement);
  h2001 = (2*ev0*eye(nphase)-jac)\(h2001-2*h20*gg1(2));
  h1101 = temp5 + hessp(:,:,2)*h11 ;
  h1101 = h1101 + multilinear3(lds.func,tens,q0,q1,h0001,x,p,ten3Increment);
  h1101 = h1101 + multilinear2(lds.func,hess,h11,h0001,x,p,hessIncrement);
  h1101 = h1101 + 2*real(multilinear2(lds.func,hess,q1,h1001,x,p,hessIncrement));
  h1101 = -jac\(h1101-2*real(gg1(2))*h11);
  h2101 = temp6 + hessp(:,:,2)*h21 ;
  h2101 = h2101 + multilinear4(lds.func,ten4,q0,q0,q1,h0001,x,p,ten4Increment);
  h2101 = h2101 + 2*multilinear3(lds.func,tens,q0,h11,h0001,x,p,ten3Increment);
  h2101 = h2101 + 2*multilinear3(lds.func,tens,q0,q1,h1001,x,p,ten3Increment);
  h2101 = h2101 + multilinear3(lds.func,tens,h20,q1,h0001,x,p,ten3Increment);
  h2101 = h2101 + multilinear3(lds.func,tens,q0,q0,conj(h1001),x,p,ten3Increment);
  h2101 = h2101 + 2*multilinear2(lds.func,hess,q0,h1101,x,p,hessIncrement);
  h2101 = h2101 + 2*multilinear2(lds.func,hess,h11,h1001,x,p,hessIncrement);
  h2101 = h2101 + multilinear2(lds.func,hess,h21,h0001,x,p,hessIncrement);
  h2101 = h2101 + multilinear2(lds.func,hess,h20,conj(h1001),x,p,hessIncrement);
  h2101 = h2101 + multilinear2(lds.func,hess,q1,h2001,x,p,hessIncrement);
  gg2= p0'*[h2110 h2101]/2;
  v01=real([gg1;gg2])\s2;
  imagp2 = imag(gg1*v01);
%generate a new point
  x0 = zeros(ntst*nphase+2,1);
  v0 = zeros(ntst*nphase+2,1);
  newbase = x+(h11-2*lyap2*(h0010*v01(1)+h0001*v01(2)))*eps^2;
for i=1:lds.ncoords/lds.nphase
    zz = exp(sqrt(-1.0)*2*pi*i/(lds.ncoords/lds.nphase -1));
    temp = newbase + real(2*q0*zz*eps + h20*zz^2*eps^2 + (h30/3*zz^3+h21*zz)*eps^3);
    x0(((i-1)*nphase+1):((i-1)*nphase+nphase)) = temp-4*lyap2*real(zz*h1010*v01(1)+zz*h1001*v01(2))*eps^3;
    temp = 2*(newbase-x)/eps+real(2*q0*zz + 2*h20*zz^2*eps+ (h30*zz^3+3*h21*zz)*eps^2);
    v0(((i-1)*nphase+1):((i-1)*nphase+nphase)) = temp-12*lyap2*real(zz*h1010*v01(1)+zz*h1001*v01(2))*eps^2;
end
% generate a new mesh and interpolate 
resc = norm(v0);
[x0,v0]=new_mesh(x0,v0,ntst,ncol);

% last entries of x0 and v0: period and parameters
omega=imag(ev0);
x0(lds.ncoords+1)     = 2*pi/(omega+(-2*lyap2*imagp2+imag(g21))*eps^2);
x0(lds.ncoords+(2:3)) = pp(ap)-2*lyap2*v01*eps^2;
v0(lds.ncoords+1)     = -4*pi*eps*(-2*lyap2*imagp2+imag(g21))/omega^2;
v0(lds.ncoords+(2:3)) = -4*lyap2*v01*eps;
v0(lds.ncoords+(1:3)) = v0(lds.ncoords+(1:3))/resc;
v0=v0/norm(v0);


% save xot x0 v0
ups = reshape(x0(lds.coords),lds.nphase,lds.tps);
[xt,p,T] = rearr(x0);
p = num2cell(p);
pars1 = lds.ncoords+1;
jac = spalloc(lds.ncoords,lds.ncoords,(lds.ncol+4)*lds.nphase);

% function
range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
for j=lds.tsts
  xp = ups(:,range0)*lds.wt;
  jac(range1,[range2 pars1]) = bordBVP_LPC_f(lds.func,xp,p,T,j);
  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end
% boundary conditions
range  = (lds.tps-1)*lds.nphase+ (lds.phases);
range1 = lds.ncoords-lds.nphase+lds.phases;
jac(range,[lds.phases range1]) = bordBVP_LPC_bc1;

range = range(lds.nphase)+1;
jac(range,lds.coords) = bordBVP_LPC_bc2(lds.func,ups,p);

%compute borders

b = []; b(lds.ncoords+2)=1; b=b';
jace=[jac,rand(lds.ncoords+1,1);rand(1,lds.ncoords+1),0];
q=jace\b;q=q(1:lds.ncoords+1,1);q=q/norm(q);
p=jace'\b;p=p(1:lds.ncoords+1,1);p=p/norm(p);

%[Q,R,E] = qr(full(jac));
%if (R(end,end) >= cds.options.FunTolerance)
%    debug('Last singular value exceeds tolerance\n');
%end
%p = E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
%lds.LPC_phi = p'/norm(p);
%p = Q(:,end);
%lds.LPC_psi = p';

lds.LPC_phi=q';
lds.LPC_psi=p';






% ---------------------------------------------------------------
function [x,p,T] = rearr(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x), parameters (p) and period (T)
global lds

nap = length(lds.ActiveParams);

p = lds.P0;
p(lds.ActiveParams) = x0(lds.ncoords+1+(1:nap));
x = x0(lds.coords);
T = x0(lds.ncoords+1);
lds.T = T;

%-----------------------------------------------------------------
function init_lds(odefile,x,p,s,ap,ntst,ncol,varargin)
global lds
lds=[];
lds.odefile = odefile;
func_handles = feval(lds.odefile);
lds.func = func_handles{2};
lds.Jacobian  = func_handles{3};
lds.JacobianP = func_handles{4};
lds.Hessians  = func_handles{5};
lds.HessiansP = func_handles{6};
lds.Der3 = func_handles{7};
lds.Der4 = func_handles{8};
lds.Der5 = func_handles{9};
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
lds.pars = lds.ncoords+(1:3);
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
if size(varargin,1)>0
    lds.BranchParams = varargin{1};
else lds.BranchParams=[];
end
lds.tsts = 1:lds.ntst;
lds.cols = 1:lds.ncol;
