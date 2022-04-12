function coef = nf_LP(odefile,jacobian,hessians,x,p,q0,p0,nphase)
%
% coef = nf_LP(odefile,jacobian,hessians,x,p,q0,p0,nphase)
% compute fold normal form coefficient.
%
global cds
  hessIncrement = (cds.options.Increment)^(3.0/4.0);
if (cds.options.SymDerivative >= 2)
  hess = chess(odefile,jacobian,hessians,x,p,ones(length(p),1));
else
  hess = [];
end
  h2 = multilinear2(odefile,hess,q0,q0,x,p,hessIncrement);	% B(q0,q0)
  coef = p0'*h2/2.0;
