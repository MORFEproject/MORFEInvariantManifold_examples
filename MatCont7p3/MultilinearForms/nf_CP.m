function coef = nf_CP(odefile,jacobian,hessians,der3,x,p,nphase)
%
% coef = nf_CP(odefile,jacobian,hessians,der3,x,p,nphase)
% compute cusp normal form coefficient.
%
global cds lpds
  jac = cjac(odefile,jacobian,x,p,lpds.ActiveParams);
  [X,D] = eig(jac);
  [val, index1] = min(abs(diag(D)));
  q0 = X(:,index1);
  [X,D] = eig(jac');
  [val, index1] = min(abs(diag(D)));
  qad0 = X(:,index1);
  p0 = qad0/(q0'*qad0);
  hessIncrement = (cds.options.Increment)^(3.0/4.0);
  ten3Increment = (cds.options.Increment)^(3.0/5.0);
if (cds.options.SymDerivative >= 3)
  hess = chess(odefile,jacobian,hessians,x,p,ones(length(p),1));
  tens = ctens3(odefile,jacobian,hessians,der3,x,p,ones(length(p),1));
else
  hess = [];
  tens = [];
end
  h2 = multilinear2(odefile,hess,q0,q0,x,p,hessIncrement);		% B(q0,q0)
  f2 = p0'*h2/2.0;
  h2 = [jac q0 ; p0' 0] \ [2*f2*q0-h2 ; 0]; 				%-A^{INV}( B(q0,q0) - <p0,B(q0,q0)>q0 ) 
  h2 = h2(1:nphase);
  h3 = multilinear3(odefile,tens,q0,q0,q0,x,p,ten3Increment);		%  C(q0,q0,q0)
  h3 = h3 +3*multilinear2(odefile,hess,q0,h2,x,p,hessIncrement);	%+3B(q0,h2)
  coef = p0'*h3/6.0;
