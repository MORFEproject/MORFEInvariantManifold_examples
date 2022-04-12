function result=lyapunov(x)
%
% calculates Lyapunov coefficient
%
global eds
[x1,p] =rearr(x);p1=num2cell(p);
jac=cjac(eds.func,eds.Jacobian,x1,p1,eds.ActiveParams);
nphase = size(x1,1);
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
if (imag(d(idx1)) == 0) && (imag(d(idx2)) == 0)
  debug('Neutral saddle\n');
  result='Neutral saddle';
  return;
end
[Q,R]=qr([real(V(:,idx1)) imag(V(:,idx1))]);
borders.v=Q(:,1:2);
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
[Q,R]=qr([real(V(:,idx1)) imag(V(:,idx1))]);
borders.w=Q(:,1:2);
k=real(d(idx1)*d(idx2));
% calculate eigenvalues
% ERROR OR WARNING
RED=jac*jac+k*eye(nphase);
Bord=[RED borders.w;borders.v' zeros(2)];
bunit=[zeros(nphase,2);eye(2)];
vext=Bord\bunit;
wext=Bord'\bunit;
result=nf_H(eds.func,eds.Jacobian,eds.Hessians,eds.Der3,x1,p1,k,vext,wext,nphase,eds.ActiveParams);

% ---------------------------------------------------------------
function [x,p] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global  eds
p = eds.P0;
p(eds.ActiveParams) = x0(end);
x = x0(1:end-1);



