function result = a_lp(x)
global eds

p = eds.P0;
p(eds.ActiveParams) = x(end);
x1 = x(1:end-1);
p1=num2cell(p);
nphase=size(x1,1);
jac=cjac(eds.func,eds.Jacobian,x1,p1,eds.ActiveParams);
% calculate eigenvalues
V=eig(full(jac));
[Y,i]=min(abs(V));
RED=jac-V(i)*eye(nphase);
% compute eigenvectors
borders.v=null(RED);
borders.w=null(RED');
Bord=[jac borders.w;borders.v' 0];
bunit=[zeros(nphase,1);1];
vext=Bord\bunit;
vext=vext(1:nphase);
wext=Bord'\bunit;
wext=wext(1:nphase);
% normalize eigenvectors
vext = vext/norm(vext);
wext = wext/(wext'*vext);
% call to normal form computation file for limit point
result = nf_LP(eds.func,eds.Jacobian,eds.Hessians,x1,p1,vext,wext,nphase);
