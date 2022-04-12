function [x0,v0]= init_BP_EP(odefile, x, p, s, h)
%
%[x0,v0]= init_BP_EP(odefile, x, p, s, h)
%
% Defines odefile
% Sets all parameters for an equilibrium continuation (p)
% and the active parameter (ap)
%
global eds cds
x0 = [x;p(eds.ActiveParams)];

eds.P0 = p;
eds.ActiveParams = eds.ActiveParams;
eds.odefile = odefile;
if ~isfield(s.data,'v') && ~(length(eds.v)==length(x0))
    [x0,v0] = init_EP_EP(odefile,x,p,eds.ActiveParams);
    return
else
    eds.v = s.data.v;
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
if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.ndim = length(x0);
cds.symjac = 1;
cds.symhess = 1;


eds.odefile = odefile;
eds.func = func_handles{2};
eds.Jacobian  = func_handles{3};
eds.JacobianP = func_handles{4};
eds.Hessians  = func_handles{5};
eds.HessiansP = func_handles{6};
eds.Der3 = func_handles{7};
eds.Der4 = func_handles{8};
eds.Der5 = func_handles{9};
siz = size(func_handles,2);
if siz > 9
    j=1;
    for i=10:siz
        eds.user{j}= func_handles{i};
        j=j+1;
    end
end
jac = cjac(cds.curve_func,cds.curve_jacobian,x0,[]);

if issparse(jac)
  Q=qr(jac);
  R=qr(jac');
  if (Q(size(x,1),size(x,1)) >= cds.options.FunTolerance)||(Q(size(x,1),size(x,1)+1) >= cds.options.FunTolerance)
    debug('Last singular value exceeds tolerance\n');
  end
  q = Q(1:end-1,1:end-2)\-Q(1:end-1,end-1:end);
  q = orth([q;eye(2)]);
  q1 = q(:,1); q2 = q(:,2);
  if (R(size(x,1),size(x,1)) >= cds.options.FunTolerance)
    debug('Last singular value exceeds tolerance\n');
  end
  p1 = [R(1:end-1,1:end-1)\-R(1:end-1,end);1];
else
  [Q,R,E] = qr(jac');
  q1 = Q(:,end-1);
  q2 = Q(:,end);
  if (R(size(x,1),size(x,1)) >= cds.options.FunTolerance)
    debug('Last singular value exceeds tolerance\n');
  end
  p1 = E*[R(1:end-1,1:end-1)\-R(1:end-1,end);1];
end

p1 = p1/norm(p1);
t1(length(x0)) = 0; t2(length(x0)) = 0;
for i=1:size(x0,1)
    x1 = x0; x1(i) = x1(i)-cds.options.Increment;
    x2 = x0; x2(i) = x2(i)+cds.options.Increment;
    t1(i) = p1'*((cjac(cds.curve_func,cds.curve_jacobian,x2,[])-cjac(cds.curve_func,cds.curve_jacobian,x1,[]))/(2*cds.options.Increment))*q1;
    t2(i) = p1'*((cjac(cds.curve_func,cds.curve_jacobian,x2,[])-cjac(cds.curve_func,cds.curve_jacobian,x1,[]))/(2*cds.options.Increment))*q2;
end
a = (t1*q1)/2;
b = (t1*q2)/2;
c = (t2*q2)/2;
if (abs(c) >= abs(a)) && (abs(c) > cds.options.FunTolerance)
    alpha1 = 1; alpha2 = 1;
    poly = [c 2*b a];
    beta = roots(poly);
    beta1 = beta(1);
    beta2 = beta(2);
elseif (abs(a) > abs(c))&&(abs(a) > cds.options.FunTolerance)
    beta1 = 1;beta2 = 1;
    poly = [a 2*b c];
    alpha = roots(poly);
    alpha1 = alpha(1);
    alpha2 = alpha(2);
else
    if(abs(b) < cds.options.FunTolerance)
        debug('Degeneracy in branching equation suspected\n');
    end
    alpha1 = 1;beta1 = 0;
    alpha2 = 0;beta2 = 1;
end
q3 = alpha1*q1+beta1*q2;
q4 = alpha2*q1+beta2*q2;
w1 = q3/norm(q3);
w2 = q4/norm(q4);
if abs(w1'*s.data.v) < abs(w2'*s.data.v)
    v0 = w1;
else
    v0 = w2;    
end
x0 = x0+h*v0;


eds.nphase = length(x);
n = eds.nphase;
a = reshape(1:(n^2),n,n);
[bia,bin,bip] = bialt(a);
if any(any(bip))
    [eds.BiAlt_M1_I,eds.BiAlt_M1_J,eds.BiAlt_M1_V] = find(bip);
else
    eds.BiAlt_M1_I=1;eds.BiAlt_M1_J=1;eds.BiAlt_M1_V=n^2+1;
end    
if any(any(bin))
    [eds.BiAlt_M2_I,eds.BiAlt_M2_J,eds.BiAlt_M2_V] = find(bin);
else
     eds.BiAlt_M2_I=1;eds.BiAlt_M2_J=1;eds.BiAlt_M2_V=n^2+1;
end
if any(any(bia))
    [eds.BiAlt_M3_I,eds.BiAlt_M3_J,eds.BiAlt_M3_V] = find(bia);
else
    eds.BiAlt_M3_I=1;eds.BiAlt_M3_J=1;eds.BiAlt_M3_V=n^2+1;
end

jac = cjac(eds.func,eds.Jacobian,x,num2cell(p),eds.ActiveParams);
jac = jac(:,1:cds.ndim-1);
jac(:,end+1) = 0;
A1 = sparse(eds.BiAlt_M1_I,eds.BiAlt_M1_J,jac(eds.BiAlt_M1_V));
A2 = sparse(eds.BiAlt_M2_I,eds.BiAlt_M2_J,jac(eds.BiAlt_M2_V));
A3 = sparse(eds.BiAlt_M3_I,eds.BiAlt_M3_J,jac(eds.BiAlt_M3_V));
A = A1-A2+A3;
[Q,R,E] = qr(full(A));
eds.bigW = Q(:,end);
eds.bigV = E(:,end);
eds.bigD = 0;
eds.adaptchoice = 1;