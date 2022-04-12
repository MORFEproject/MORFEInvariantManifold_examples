function [x,v] = init_BT_Hom(odefile, x, s, p, ap, ntst, ncol,TTolerance ,amplitude, extravec)
global homds cds
%
% 1. CHECCK INPUT
%--------------------------------------------
n_par = size(ap,2);
if n_par ~= 2
    error('2 free system parameters are needed');
end
%
if extravec==zeros(1,3)
    error('at least 1 free homoclinic parameter is needed');
end
%
if extravec==ones(1,3)
    error('at most 2 free homoclinic parameters are needed');
end
%
for i=1:length(s)
    if ~isempty(s.data) && isfield(s(i).data,'evec')
        nph = size(s(i).data.evec,2);
    elseif isfield(s(i).data,'x')
        nph = size(s(i).data.x,2);
    else
        
    end
end
%
if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end
cds.curve = @homoclinic;
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_jacobian = curvehandles{4};
cds.curve_hessians = curvehandles{5};
homds = [];
%
% 2. INITIALIZE homds
%--------------------------------------------
eps0=1e-2;
eps1=1e-2;
   T=100;
%
init_homds(odefile,x,p,ap,ntst,ncol,extravec,T,eps0,eps1,nph);
%
func_handles = feval(odefile);
symord = 0; 
symordp = 0;
%
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
%
homds.odefile = odefile;
homds.func = func_handles{2};
homds.Jacobian  = func_handles{3};
homds.JacobianP = func_handles{4};
homds.Hessians  = func_handles{5};
homds.HessiansP = func_handles{6};
homds.Der3 = func_handles{7};
homds.Der4 = func_handles{8};
homds.Der5 = func_handles{9};
%
cds.oldJac = [];
cds.oldJacX = [];
xp = [x;p(ap)];
cds.ndim = length(xp);
%
homds.x0 = x(1:homds.nphase);
%  
pcell = num2cell(p);
%
% 3. COMPUTE A, J1
%--------------------------------------------
A = cjac(homds.func,homds.Jacobian,homds.x0,pcell,homds.ActiveParams);
J1 = cjacp(homds.func,homds.JacobianP,homds.x0,pcell,ap);
%
% 4. COMPUTE B, A1, J2, q0, q1, p0 and p1
%--------------------------------------------
% i)  B & A1
B =chess(homds.func,homds.Jacobian,homds.Hessians,homds.x0,pcell,homds.ActiveParams);
A1=chessp(homds.func,homds.Jacobian,homds.HessiansP,homds.x0,pcell,homds.ActiveParams);
if size(A1,3) > length(ap)
    A1 = A1(:,:,ap);
end
%
% ii) J2
for i=homds.ActiveParams
    pa1 = pcell; pa1{i} = pa1{i}-cds.options.Increment;
    pa2 = pcell; pa2{i} = pa2{i}+cds.options.Increment;
    Hjp2 = cjacp(homds.func,homds.JacobianP,homds.x0,pa2,ap); 
    Hjp1 = cjacp(homds.func,homds.JacobianP,homds.x0,pa1,ap); 
    tmpJ2(:,:,i) = Hjp2 - Hjp1;
end
tmpJ2 = tmpJ2/(2*cds.options.Increment);
J2 = tmpJ2;
if size(J2,3) > length(ap)
    J2 = J2(:,:,ap);
end
% iii) q0, q1, p0 & p1
  [X,D] = eig(A);
  index1 = find(abs(diag(D)) < 1e-3); %If ok, index1 is 1x2 array otherwise
  vext = real(X(:,index1(1)));
  [X,D] = eig(A');
  index1 = find(abs(diag(D)) < 1e-3);
  wext = real(X(:,index1(1)));
  Bord = [ A wext; vext' 0];
  bunit=[zeros(homds.nphase,1);1];
  q0=Bord\bunit; 
  q0=q0(1:homds.nphase);          % A q0 = 0 , <vext,q0> = 1
  p1=Bord'\bunit;
  p1=p1(1:homds.nphase);          % A'p1 = 0 , <wext,p1> = 1
  Bord = [ A p1; q0' 0];
  q1 = Bord\[q0; 0];		
  q1 = q1(1:homds.nphase);		  % A q1 = q0, <q0, q1>  = 0
  p0 = Bord'\[p1; 0];
  p0 = p0(1:homds.nphase);		  % A'p0 = p1, <p0, p1>  = 0
%
% 5. NORMALIZE SO THAT <p0,q0>=<p1,q1>=1 & <p0,q1>=<p1,q0>=0 
%--------------------------------------------
  mu = sqrt(abs(q0'*q0));
  q0 = (1/mu)*q0;
  q1 = (1/mu)*q1;
  q1 = q1 - (q0'*q1)*q0;
  nu = q0'*p0;
  p1 = (1/nu)*p1;
  p0 = p0 - (p0'*q1)*p1;
  p0 = (1/nu)*p0;
%
% 6. BT NORMAL FORM COEFFICIENTS (a,b)
%--------------------------------------------
hessIncrement=(cds.options.Increment)^(3.0/4.0);
ten3Increment=(cds.options.Increment)^(3.0/5.0);

if (cds.options.SymDerivative>=3)
  hess=chess(homds.func,homds.Jacobian,homds.Hessians,homds.x0,pcell,homds.ActiveParams);
  tens=ctens3(homds.func,homds.Jacobian,homds.Hessians,homds.Der3,homds.x0,pcell,homds.ActiveParams);
else
  hess = [];
  tens = [];
end
M1 = multilinear2(homds.func,hess,q0,q0,homds.x0,pcell,hessIncrement);	% B(q0,q0)
M2 = multilinear2(homds.func,hess,q0,q1,homds.x0,pcell,hessIncrement);	% B(q0,q1)
M3 = multilinear2(homds.func,hess,q1,q0,homds.x0,pcell,hessIncrement);	% B(q1,q0)
M4 = multilinear2(homds.func,hess,q1,q1,homds.x0,pcell,hessIncrement);	% B(q1,q1)
% 
a = 0.5*p1'* M1;
b = p1'*M2+p0'*M1;
%
% 7.(a) COMPUTE H2000, H1100, H0200, H00, K1, K22, H1001, H0002 and H0101 (THE HOMOCLINIC PREDICTOR COEFFICIENTS)
%----------------------------------------------------------------------------------------------------------------
 for i=1:homds.nphase
    M5(:,i)=p1'*B(:,:,i)*q1;
    M6(:,i)=p1'*B(:,:,i)*q0;
    M7(:,i)=p0'*B(:,:,i)*q0;
    M8(:,i)=p0'*B(:,:,i)*q1;    
end
for i=1:length(ap)
     M9(:,i)=p1'*A1(:,:,i)*q1;
    M10(:,i)=p1'*A1(:,:,i)*q0;
    M11(:,i)=p0'*A1(:,:,i)*q0;
    M12(:,i)=p0'*A1(:,:,i)*q1;  
end
%
H2000=Bord\[2*a*q1-M1; 0];
H2000=H2000(1:end-1,:);
%
% Correction by Maikel Bosschaert to enure that H0200 is solvable
H2000=H2000+1/2*(-2*p0'*H2000 + 2*p0'*M2 + p1'*M4)*q0;
%
H1100=Bord\[b*q1-M2+H2000; 0];
H1100=H1100(1:end-1,:);
%
H0200=Bord\[2*H1100-M4; 0];
H0200=H0200(1:end-1,:);
%
    C=[M6; M7+M5];             D=[M10; M11+M9];
    F=[q1 zeros(size(A,1),1)]; G=[0.5*p1'*M4 0; -p0'*M4+3*p0'*H1100 1];
tmpHK=[A J1; C D];
   HK=tmpHK\[F; G];
%
H00=HK(1:end-2,:);   
 K1=HK(end-1:end,:);
%
for i = 1:homds.nphase
    N1(:,i)=p1'*B(:,:,i)*H00(:,2);
end
for i=1:length(ap)
    N2(:,i)=p1'*A1(:,:,i)*H00(:,2);
    N3(:,i)=p1'*J2(:,:,i)*K1(:,2) ;
end
N1=N1* H00(:,2);
N2=N2*K1(:,2);
N3=N3*K1(:,2);
K22=-(N1+2*N2+N3)*K1(:,1);
%--------------------------------------------------------------------------
M13=multilinear2(homds.func,hess,q0,H00(:,2),homds.x0,pcell,hessIncrement);	       % B(q0,H0001)
%
for i=1:length(ap)
    M14(:,i)=A1(:,:,i)*q0;
    M15(:,i)=A1(:,:,i)*q1;
end
%
M14=M14*K1(:,2);   % A1(q0,K1.1)
M15=M15*K1(:,2);   % A1(q1,K1.1)
%
M16=multilinear2(homds.func,hess,H00(:,2),H00(:,2),homds.x0,pcell,hessIncrement);  % B(H0001,H0001)
%
for i=1:length(ap)
    M17(:,i)=A1(:,:,i)*H00(:,2); 
    M18(:,i)=J2(:,:,i)*K1(:,2);
end
M17=M17*K1(:,2);   %A1(H0001,K1.1)
M18=M18*K1(:,2);   %J2(K1.1,K1.1)
  z=M16+2*M17+M18;
%
H1001=Bord\[-(M13+M14) ; 0];
H1001=H1001(1:end-1,:);
%
H0002=Bord\[-(z+J1*K22); 0];
H0002=H0002(1:end-1,:);
%
M19=multilinear2(homds.func,hess,q1,H00(:,2),homds.x0,pcell,hessIncrement);	       % B(q1,H0001)
%
H0101=Bord\[-(M19+M15-H1001-q1); 0];
H0101=H0101(1:end-1,:);
%
% 7.(b): COMPUTE d, e, a1 AND b1 (THE BT NORMAL FORM COEFFICIENTS)
%-----------------------------------------------------------------
%
for i=homds.ActiveParams
    pa1= pcell; pa1{i} = pa1{i}-cds.options.Increment;
    pa2= pcell; pa2{i} = pa2{i}+cds.options.Increment;
    Hp2=chess(homds.func,homds.Jacobian,homds.Hessians,homds.x0,pa2,homds.ActiveParams);
    Hp1=chess(homds.func,homds.Jacobian,homds.Hessians,homds.x0,pa1,homds.ActiveParams);
    %
    B3(:,:,:,i)=Hp2-Hp1; 
end
B3=B3/(2*cds.options.Increment);
%
if size(B3,4) > length(ap)
   B3=B3(:,:,:,ap);
end
%
for i=1:length(ap)
    B2=B3(:,:,:,i);
    tmpM20=tensor2op(B2,q0,q0,homds.nphase);
    tmpM21=tensor2op(B2,q0,q1,homds.nphase);
    M20(:,i)=tmpM20;
    M21(:,i)=tmpM21;
end
M20=M20*K1(:,2);   %B2(q0,q0,K1.1);
M21=M21*K1(:,2);   %B2(q0,q1,K1.1);
%
M22=multilinear3(homds.func,tens,q0,q0,q0,homds.x0,pcell,ten3Increment);        %C(q0,q0,q0)
M23=multilinear3(homds.func,tens,q0,q0,q1,homds.x0,pcell,ten3Increment);        %C(q0,q0,q1)
M24=multilinear3(homds.func,tens,q0,q1,H00(:,2),homds.x0,pcell,ten3Increment);  %C(q0,q1,H0001)
M25=multilinear3(homds.func,tens,q0,q0,H00(:,2),homds.x0,pcell,ten3Increment);  %C(q0,q0,H0001)
%
M26=multilinear2(homds.func,hess,q0,H2000,homds.x0,pcell,hessIncrement);        %B(q0,H2000)
M27=multilinear2(homds.func,hess,q0,H1100,homds.x0,pcell,hessIncrement);        %B(q0,H1100)
M28=multilinear2(homds.func,hess,q1,H2000,homds.x0,pcell,hessIncrement);        %B(q1,H2000)
M29=multilinear2(homds.func,hess,q0,H1001,homds.x0,pcell,hessIncrement);        %B(q0,H1001)
M30=multilinear2(homds.func,hess,H00(:,2),H2000,homds.x0,pcell,hessIncrement);  %B(H0001,H2000)
M31=multilinear2(homds.func,hess,q1,H1001,homds.x0,pcell,hessIncrement);        %B(q1,H1001)
M32=multilinear2(homds.func,hess,H00(:,2),H1100,homds.x0,pcell,hessIncrement);  %B(H0001,H1100)
M33=multilinear2(homds.func,hess,q0,H0101,homds.x0,pcell,hessIncrement);        %B(q0,H0101)
%
for i=1:length(ap)
   M34(:,i)=A1(:,:,i)*H2000;
   M35(:,i)=A1(:,:,i)*H1100;
end
M34=M34*K1(:,2);   %A1(H2000,K1.1)
M35=M35*K1(:,2);   %A1(H1100,K1.1)
%
d=p1'*(1/6*M22+0.5*M26-a*H1100);
%
H3000=Bord\[-6*(1/6*M22+0.5*M26-a*H1100-d*q1); 0];
H3000=H3000(1:end-1,:);
%
e=p1'*(0.5*M23+M27+0.5*M28-b*H1100-a*H0200-0.5*H3000);
%
a1=p1'*(0.5*M25+0.5*M20+M29+0.5*M30+0.5*M34-a*H0101);
%
H2001=Bord\[-2*(0.5*M25+0.5*M20+M29+0.5*M30+0.5*M34-a*H0101-a1*q1); 0];
H2001=H2001(1:end-1,:);
%
b1=p1'*(M24+M21+M31+M32+M33+M35-b*H0101-H1100-H2001);
%
eps=sqrt(amplitude*abs(a)/6);
%



%
% 8. INITIAL CYCLE
%--------------------------------------------
% a) THE INITIAL  APPROXIMATION OF ALPHA  
%
alpha=((10*b)/(7*a))*eps^2*K1(:,2)+eps^4*((-4/a)*K1(:,1)+((50*b^2)/(49*a^2))*K22+...
    b/a*K1(:,2)*(1/a*(100/49*b1-4*e/b)+1/a^2*(-50/49*b*a1+288/2401*b^2+146/49*d)));
%
p(ap) = p(ap) + alpha;  % tmpfreep=tmpfreep+alpha;
homds.P0 = p;
%
% b) THE INITIAL HALF-RETURN TIME VALUE T:-
% |w0(-inf,eps)-w0(T,eps)|=delta0=>(6*norm(q0)*eps^2/|a|)*tanh(eps*T)^2=k
% where norm(q0)=1.
myarg = sqrt((TTolerance * abs(a)) / (6 * eps^2));
myargn = log(1/myarg + ((1/myarg)^2 - 1));
if ~isempty(imag(myargn)) && (imag(myargn) ~= 0)
    myargn = log(1/myarg - ((1/myarg)^2 - 1));
end
%
homds.T = abs(myargn/eps);
%

disp('BT normal form coefficients:')
%disp(table(a,b,d,e,a1,b1));
disp(table(a,b));  %NN: removed d,e,a1,b1 as requested.

fprintf('T: %g\n', homds.T);
name = 'The initial perturbation parameter:';
fprintf('%s  %d\n',name,eps);


% c) THE INITIAL APPROXIMATION OF CYCLE
%
homds.x0 = homds.x0 + eps^2 * ((10*b)/(7*a) * H00(:,2) + 2/a * q0);
x0 = homds.x0;
for i=1:length(homds.finemsh)    
    % CONVERSION FROM [0,1] --> [-T,+T]
    t = (2*homds.finemsh(i) - 1) * (homds.T);
    x0 = homds.x0;
    %
    ups(:,i)=((10/7)*H00(:,2)*b/a + q0*(-6*sech(eps*t)^2+2)/a)*eps^2+...
        12*q1*sech(eps*t)^2*tanh(eps*t)*eps^3/a + (q0*(-(1/49)*(-210*...
        a1*b+18*b^2+ 147*d)*sech(eps*t)^2/a^2-(2/7)*(5*a1*b+7*d)/a^2)...
        /a- (72/7)*q1*b*sinh(eps*t)^2/(a^2*cosh(eps*t)^4)+(10/7)*H1001...
        *(-6*sech(eps*t)^2+2)*b/a^2 +(50/49)*H0002*b^2/a^2+(1/2)*H2000*...
        (-6*sech(eps*t)^2+2)^2/a^2+H00(:,2)*((100/49)*b1/a-(50/49)*a1...
        *b/a^2+(288/2401)*b^2/a^2- 4*e/(a*b) + (146/49)*d/a^2)*b/a-4*...
        H00(:,1)/a)*eps^4; 
    %
    ups(:,i) = ups(:,i) + homds.x0;
end
%
%plot(ups(1,:),ups(2,:))
%
%d) THE INITIAL eps0 & eps1 :
%
homds.eps0 = norm(ups(:,1) - homds.x0);
homds.eps1 = norm(ups(:,end) - homds.x0);
%
% 9. DIMENSION OF STABLE & UNSTABLE SUBSPACE 
%--------------------------------------------
Hom_calc_weights;
A = cjac(homds.func,homds.Jacobian,x0,num2cell(p),homds.ActiveParams);
D = eig(A);
homds.nneg = sum(real(D) < 0);
homds.npos = homds.nphase-homds.nneg;
homds.Ysize = homds.nneg*homds.npos;
%
% 10. COMPOSE x
%--------------------------------------------
% 
% i. CYCLE
%
x1 = reshape(ups,size(ups,1)*size(ups,2),1);
v = []; 
[x1,v]=Hom_new_mesh(x1,v,ntst,ncol);
%
% ii. EQUILIBRIUM COORDINATES
%
x1 = [x1; x0];
%
% iii. TWO FREE PARAMETERS
%
x1 = [x1; homds.P0(homds.ActiveParams)];
%
% iv. EXTRA FREE PARAMETER (FREE HOMOCLINIC PARAMETER)
%
extravec = [homds.T; homds.eps0; homds.eps1];
x1 = [x1; extravec(find(homds.extravec))];

for i=1:homds.nneg
    x1 = [x1; zeros(homds.npos,1)];
end

for i=1:homds.npos
    x1 = [x1; zeros(homds.nneg,1)];
end
%
%
% 11. ASSIGN SOME VALUES TO HOMOCLINIC FIELDS
%--------------------------------------------
% a) YS AND YU, INITIALIZED TO 0
%
homds.YS = zeros(homds.npos,homds.nneg);
homds.YU = zeros(homds.nneg,homds.npos);

% b) THIRD PARAMETER = UNSTABLE_FLAG:-
%
% i) 1 IF WE WANT THE UNSTABLE SPACE. 
% ii) 0 IF WE WANT THE STABLE ONE.
[QU, se] = computeBase(A,0,homds.npos);
[QS, se] = computeBase(A,1,homds.nneg);

homds.oldStableQ = QS;
homds.oldUnstableQ = QU;
homds.ups = [];
homds.ndim = length(x1);
cd.ndim = homds.ndim;
x = x1;
%
%
% 12. INITIAL TANGENT VECTOR
%--------------------------------------------
ups0 = reshape(x(1:size(ups,1)*size(ups,2),:),homds.nphase,homds.tps);
pp1 = num2cell(p);
 for i=1:homds.tps
    homds.upoldp(:,i) = 2*T*feval(homds.func, 0, ups0(:,i), pp1{:});
 end
BTHomjac = BVP_Hom_jac(homds.func,x(1:size(ups,1)*size(ups,2),1),x0,p,homds.T,homds.eps0,homds.eps1,homds.YS,homds.YU);
[Q,R] = qr(full(BTHomjac)');
v = full(Q(:,end));

%---------------------------------------------------------------------------
function init_homds(odefile,x,p,ap,ntst,ncol,extravec,T,eps0,eps1,nph)
global homds 

homds.odefile = odefile;
func_handles = feval(homds.odefile);
homds.func = func_handles{2};
homds.Jacobian  = func_handles{3};
homds.JacobianP = func_handles{4};
homds.Hessians  = func_handles{5};
homds.HessiansP = func_handles{6};
homds.Der3=[];
siz = size(func_handles,2);

if siz > 9
    j=1;
    for k=10:siz
        homds.user{j}= func_handles{k};
        j=j+1;
    end
else homds.user=[];
end

homds.nphase = nph;
homds.ActiveParams = ap;
homds.P0 = p;
homds.extravec = extravec;

Hom_set_ntst_ncol(ntst,ncol,(0:ntst)/ntst);

homds.T = T;
homds.eps0 = eps0;
homds.eps1 = eps1;
homds.cols_p1 = 1:(homds.ncol+1);
homds.cols_p1_coords = 1:(homds.ncol+1)*homds.nphase;
homds.ncol_coord = homds.ncol*homds.nphase;
homds.col_coords = 1:homds.ncol*homds.nphase;
homds.pars = homds.ncoords+(1:3);
homds.phases = 1:homds.nphase;
homds.ntstcol = homds.ntst*homds.ncol;
homds.wp = kron(homds.wpvec',eye(homds.nphase));
homds.pwwt = kron(homds.wt',eye(homds.nphase));
homds.pwi = homds.wi(ones(1,homds.nphase),:);

homds.bialt_M1 = [];
homds.bialt_M2 = [];
homds.bialt_M3 = [];
homds.bialt_M4 = [];
homds.multipliers = nan;
homds.monodromy = [];
homds.multi_r1 = [];
homds.multi_r2 = [];
homds.ups = [];
homds.vps = [];
homds.tsts = 1:homds.ntst;
homds.cols = 1:homds.ncol;

homds.HTPstep = 0;