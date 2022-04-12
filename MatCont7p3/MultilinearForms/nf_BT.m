function coef = nf_BT(odefile,jacobian,hessians,x,p,nphase)
%
% coef = nf_bt(odefile,jacobian,hessians,x,p,vext,wext,nphase)
% compute normal form coefficients for Bogdanov-Takens.
global cds
  jac=cjac(odefile,jacobian,x,p,ones(length(p),1));
  [X,D] = eig(jac);
  index1 = find(abs(diag(D)) < 1e-3);%If ok, index1 is 1x2 array otherwise
if(isempty(index1)==1)
  debug('NF-computation impossible\n');
  coef =[ 0,0 ];
  return;
end
% MATLAB will fail to find to do a proper [X,D] = eig(jac) decomposition,
% therefore it is necessarty to use vext/wext to find (generalized) eigenvectors.
  vext = real(X(:,index1(1)));
  [X,D] = eig(jac');
  index1 = find(abs(diag(D)) < 1e-3);
  wext = real(X(:,index1(1)));
  Bord = [ jac wext; vext' 0];
  bunit=[zeros(nphase,1);1];
  q0=Bord\bunit; 
  q0=q0(1:nphase);          % A q0 = 0, <vext,q0> = 1
  p1=Bord'\bunit;
  p1=p1(1:nphase);          % A'p1 = 0, <wext,p1> = 1
  Bord = [ jac p1; q0' 0];
  q1 = Bord\[q0; 0];		
  q1 = q1(1:nphase);		% A q1 = q0, <q0,q1> = 0
  p0 = Bord'\[p1; 0];
  p0 = p0(1:nphase);		% A'p0 = p1, <p0,p1> = 0

% normalize so that <p0,q0>=<p1,q1>=1, <p0,q1>=<p1,q0>=0 and <q0,q0>=1, <q0,q1>=0
  mu = sqrt(abs(q0'*q0));
  q0 = (1/mu)*q0;
  q1 = (1/mu)*q1;
  q1 = q1 - (q0'*q1)*q0;
  nu = q0'*p0;
  p1 = (1/nu)*p1;
  p0 = p0 - (p0'*q1)*p1;
  p0 = (1/nu)*p0;

  hessIncrement = (cds.options.Increment)^(3.0/4.0);
if (cds.options.SymDerivative >= 2)
  hess = chess(odefile,jacobian,hessians,x,p,ones(length(p),1));
else
  hess = [];
end
  h20 = multilinear2(odefile,hess,q0,q0,x,p,hessIncrement);	% B(q0,q0)
  h11 = multilinear2(odefile,hess,q0,q1,x,p,hessIncrement);	% B(q0,q1)
  coef = [ p1'*h20/2.0 (p0'*h20 + p1'*h11) ];
