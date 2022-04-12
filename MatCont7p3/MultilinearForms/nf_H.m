function coef = nf_H(odefile,jacobian,hessians,der3,x,p,k,vext,wext,nphase,ap)
%
%calculates the first lyapunov coefficient
% 
global cds
%compute jacobian and eigenvectors.
  jac=cjac(odefile,jacobian,x,p,ap);
  omega=sqrt(k);
  alpha=vext(1:nphase,1)'*jac*vext(1:nphase,2)-1i*omega*vext(1:nphase,1)'*vext(1:nphase,2);
  beta=-vext(1:nphase,1)'*jac*vext(1:nphase,1)+1i*omega*vext(1:nphase,1)'*vext(1:nphase,1);
  q0=alpha*vext(1:nphase,1)+beta*vext(1:nphase,2);
  alpha=wext(1:nphase,1)'*jac'*wext(1:nphase,2)+1i*omega*wext(1:nphase,1)'*wext(1:nphase,2);
  beta=-wext(1:nphase,1)'*jac'*wext(1:nphase,1)-1i*omega*wext(1:nphase,1)'*wext(1:nphase,1);
  p0=alpha*wext(1:nphase,1)+beta*wext(1:nphase,2);
  q0=q0/norm(q0);
  p0=p0/(q0'*p0);
  hessIncrement = (cds.options.Increment)^(3.0/4.0);
  ten3Increment = (cds.options.Increment)^(3.0/5.0);
if (cds.options.SymDerivative >= 3)
  hess = chess(odefile,jacobian,hessians,x,p,ones(length(p),1));
  tens = ctens3(odefile,jacobian,hessians,der3,x,p,ones(length(p),1));
else
  hess = [];
  tens = [];
end
%2nd order vectors
  h20 = (2*i*omega*eye(nphase)-jac)\multilinear2(odefile,hess,q0,q0,x,p,hessIncrement);	% (2iw-A)\B(q0,q0)
  h11 = -jac\multilinear2(odefile,hess,q0,conj(q0),x,p,hessIncrement);			% -A\B(q0,conj(q0))
%3rd order vectors
  h21 = multilinear3(odefile,tens,q0,q0,conj(q0),x,p,ten3Increment);			%  C(q0,q0,conj(q0))
  h21 = h21 + 2*multilinear2(odefile,hess,q0,h11,x,p,hessIncrement);			%+2B(h11,q0)
  h21 = h21 + multilinear2(odefile,hess,h20,conj(q0),x,p,hessIncrement);		%+ B(h20,conj(q0))
  coef = real(p0'*h21)/2.0;
