function [x0,v0]= init_H_H(odefile, x, p, ap)
%
% [x1,v1] = init_H_H(odefile, x, p, ap)
%
% Initializes a Hopf bifurcation continuation from a Hopf point
% 

global cds hds

% check input
if size(ap,2)~=2
    errordlg('Two active parameter are needed for a Hopfpoint bifurcation continuation');
end
v0=[];
% initialize hds
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


hds.odefile = odefile;
hds.func = func_handles{2};
hds.Jacobian  = func_handles{3};
hds.JacobianP = func_handles{4};
hds.Hessians  = func_handles{5};
hds.HessiansP = func_handles{6};
hds.Der3 = func_handles{7};
hds.Der4 = func_handles{8};
hds.Der5 = func_handles{9};
siz = size(func_handles,2);
if siz > 9
    j=1;
    for i=10:siz
        hds.user{j}= func_handles{i};
        j=j+1;
    end
end
hds.nphase = size(x,1);
hds.ActiveParams = ap;
hds.P0 = p;
cds.curve = @hopf;
cds.ndim = length(x)+3;
cds.symjac = 1;
cds.symhess = 0;
if size(hds.P0,2) == 1
    x0=[x;hds.P0(ap)];
else
    x0=[x;hds.P0(ap)'];
end
[x1,p] = rearr(x0); p = num2cell(p);
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_options = curvehandles{3};
cds.curve_jacobian = curvehandles{4};
cds.curve_hessians = curvehandles{5};
cds.options = feval(cds.curve_options);
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.options = contset(cds.options,'Increment',1e-5);

jac = cjac(hds.func,hds.Jacobian,x1,p,hds.ActiveParams);
nphase = size(x1,1);
nap = length(hds.ActiveParams);
% calculate eigenvalues and eigenvectors

[V,D] = eig(jac);
% find pair of complex eigenvalues
d = diag(D);
smallest_sum = Inf;
for j=1:nphase-1
  [val,idx] = min(abs(d(j+1:nphase)+d(j)));
  if val < smallest_sum
    idx1 = j;
    idx2 = j+idx;
    smallest_sum = val;
  end
end
% real part? Oh dear, a neutral saddle!
if (imag(d(idx1)) == 0) && (imag(d(idx2)) == 0)
    [Q,R,E] = qr([V(:,idx1) V(:,idx2)]);
else
    [Q,R,E] = qr([real(V(:,idx1)) imag(V(:,idx1))]);
end
hds.borders.v = Q(:,1:2);

[V,D] = eig(jac');
% find pair of complex eigenvalues
d = diag(D);
smallest_sum = Inf;
for j=1:nphase-1
  [val,idx] = min(abs(d(j+1:nphase)+d(j)));
  if val < smallest_sum
    idx1 = j;
    idx2 = j+idx;
    smallest_sum = val;
  end
end
if (imag(d(idx1)) == 0) && (imag(d(idx2)) == 0)
    [Q,R,E] = qr([V(:,idx1) V(:,idx2)]);
else
    [Q,R,E] = qr([real(V(:,idx1)) imag(V(:,idx1))]);
end
hds.borders.w = Q(:,1:2);

k  = real(d(idx1)*d(idx2));
x0 = [x0;k];

% calculate eigenvalues
% ERROR OR WARNING
RED  = jac*jac+k*eye(hds.nphase);
A = [jac  cjacp(hds.func,hds.JacobianP,x1,p,hds.ActiveParams) zeros(hds.nphase,1)];
[Q,R] = qr(A');
Bord  = [RED hds.borders.w;hds.borders.v' zeros(2)];
bunit = [zeros(hds.nphase,2);eye(2)];
vext  = Bord\bunit;
wext  = Bord'\bunit;
hess  = chess(hds.func,hds.Jacobian,hds.Hessians,x1,p,hds.ActiveParams);
gx(4,hds.nphase) = 0;
for i = 1:hds.nphase
    gx(1,i) = -wext(1:hds.nphase,1)'*(jac*hess(:,:,i)+hess(:,:,i)*jac)*vext(1:hds.nphase,1);
    gx(2,i) = -wext(1:hds.nphase,1)'*(jac*hess(:,:,i)+hess(:,:,i)*jac)*vext(1:hds.nphase,2);
    gx(3,i) = -wext(1:hds.nphase,2)'*(jac*hess(:,:,i)+hess(:,:,i)*jac)*vext(1:hds.nphase,1);
    gx(4,i) = -wext(1:hds.nphase,2)'*(jac*hess(:,:,i)+hess(:,:,i)*jac)*vext(1:hds.nphase,2);
end
gk(1,1) = -wext(1:hds.nphase,1)'*vext(1:hds.nphase,1);
gk(2,1) = -wext(1:hds.nphase,1)'*vext(1:hds.nphase,2);
gk(3,1) = -wext(1:hds.nphase,2)'*vext(1:hds.nphase,1);
gk(4,1) = -wext(1:hds.nphase,2)'*vext(1:hds.nphase,2);
hessp = chessp(hds.func,hds.Jacobian,hds.HessiansP,x1,p,hds.ActiveParams);
gp(4,nap) = 0;
for i = 1:nap
    gp(1,i) = -wext(1:hds.nphase,1)'*(jac*hessp(:,:,i)+hessp(:,:,i)*jac)*vext(1:hds.nphase,1);
    gp(2,i) = -wext(1:hds.nphase,1)'*(jac*hessp(:,:,i)+hessp(:,:,i)*jac)*vext(1:hds.nphase,2);
    gp(3,i) = -wext(1:hds.nphase,2)'*(jac*hessp(:,:,i)+hessp(:,:,i)*jac)*vext(1:hds.nphase,1);
    gp(4,i) = -wext(1:hds.nphase,2)'*(jac*hessp(:,:,i)+hessp(:,:,i)*jac)*vext(1:hds.nphase,2);
end
A = [A;gx gp gk]*Q;
Jres = A(1+hds.nphase:end,1+hds.nphase:end)';
[Q,R,E] = qr(Jres);
index = [1 1;1 2;2 1;2 2];
[I,J] = find(E(:,1:2));
hds.index1 = index(I(1),:);
hds.index2 = index(I(2),:);
%cds = rmfield(cds,'options');

if hds.nphase > 3
    jac(:,end+1) = 0;
    A1 = sparse(hds.BiAlt_M1_I,hds.BiAlt_M1_J,jac(hds.BiAlt_M1_V));
    A2 = sparse(hds.BiAlt_M2_I,hds.BiAlt_M2_J,jac(hds.BiAlt_M2_V));
    A3 = sparse(hds.BiAlt_M3_I,hds.BiAlt_M3_J,jac(hds.BiAlt_M3_V));
    A = A1-A2+A3;
    [Q,R,E] = qr(full(A));
    hds.bigW = [Q(:,end) Q(:,end-1)];
    hds.bigV = [E(:,end) E(:,end-1)];
    hds.bigD = zeros(2,2);
    hds.adaptchoice = 1;
end

% ---------------------------------------------------------------
function [x,p] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global hds
p = hds.P0;
p(hds.ActiveParams) = x0((hds.nphase+1):end);
x = x0(1:hds.nphase);

n = hds.nphase;
a = reshape(1:(n^2),n,n);
[bia,bin,bip] = bialt(a);
if any(any(bip))
    [hds.BiAlt_M1_I,hds.BiAlt_M1_J,hds.BiAlt_M1_V] = find(bip);
else
    hds.BiAlt_M1_I=1;hds.BiAlt_M1_J=1;hds.BiAlt_M1_V=n^2+1;
end    
if any(any(bin))
    [hds.BiAlt_M2_I,hds.BiAlt_M2_J,hds.BiAlt_M2_V] = find(bin);
else
     hds.BiAlt_M2_I=1;hds.BiAlt_M2_J=1;hds.BiAlt_M2_V=n^2+1;
end
if any(any(bia))
    [hds.BiAlt_M3_I,hds.BiAlt_M3_J,hds.BiAlt_M3_V] = find(bia);
else
    hds.BiAlt_M3_I=1;hds.BiAlt_M3_J=1;hds.BiAlt_M3_V=n^2+1;
end
