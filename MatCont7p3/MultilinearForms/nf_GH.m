function coef = nf_GH(odefile,jacobian,hessians,der3,der4,der5,x,p,k,vext,wext,nphase)

%
%calculates the second lyapunov coefficient
% 
global cds hds
%compute jacobian and eigenvectors.
  jac=cjac(odefile,jacobian,x,p,hds.ActiveParams);
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
  ten4Increment = (cds.options.Increment)^(3.0/6.0);
  ten5Increment = (cds.options.Increment)^(3.0/7.0);
if (cds.options.SymDerivative >= 3)
  hess = chess(odefile,jacobian,hessians,x,p,ones(length(p),1));
  tens = ctens3(odefile,jacobian,hessians,der3,x,p,ones(length(p),1));
else
  hess = [];
  tens = [];
end
if (cds.options.SymDerivative >= 5)
  ten4 = feval(hds.Der4,0,x,p{:});
  ten5 = feval(hds.Der5,0,x,p{:});
else
  ten4 = [];
  ten5 = [];
end

%2nd order vectors
  h20 = (2i*omega*eye(nphase)-jac)\multilinear2(odefile,hess,q0,q0,x,p,hessIncrement);	% (2iw-A)\B(q0,q0)
  h11 = -jac\multilinear2(odefile,hess,q0,conj(q0),x,p,hessIncrement);			% -A\B(q0,conj(q0))
%3rd order vectors
  h30 = multilinear3(odefile,tens,q0,q0,q0,x,p,ten3Increment);				%  C(q0,q0,conj(q0))
  h30 = h30 + 3*multilinear2(odefile,hess,q0,h20,x,p,hessIncrement);			%+3B(h20,q0)  
  h30 = (3i*omega*eye(nphase)-jac)\h30;
  h21 = multilinear3(odefile,tens,q0,q0,conj(q0),x,p,ten3Increment);			%  C(q0,q0,conj(q0))
  h21 = h21 + 2*multilinear2(odefile,hess,q0,h11,x,p,hessIncrement);			%+2B(h11,q0)
  h21 = h21 + multilinear2(odefile,hess,h20,conj(q0),x,p,hessIncrement);		%+ B(h20,conj(q0))
  g21 = p0'*h21/2.0;
  h21 = [jac-1i*omega*eye(nphase) q0; p0' 0 ]\[ 2*g21*q0-h21 ; 0];
  h21 = h21(1:nphase);
%4th order vectors
  h31 = multilinear4(odefile,ten4,q0,q0,q0,conj(q0),x,p,ten4Increment);			%  D(q0,q0,q0,conj(q0))
  h31 = h31 + 3*multilinear3(odefile,tens,q0,q0,h11,x,p,ten3Increment);			%+3C(q0,q0,h11)
  h31 = h31 + 3*multilinear3(odefile,tens,q0,conj(q0),h20,x,p,ten3Increment);		%+3C(q0,conj(q0),h20)
  h31 = h31 + 3*multilinear2(odefile,hess,h20,h11,x,p,hessIncrement);			%+3B(h20,h11)
  h31 = h31 + 3*multilinear2(odefile,hess,h21,q0,x,p,hessIncrement);			%+3B(h21,q0)
  h31 = h31 +   multilinear2(odefile,hess,h30,conj(q0),x,p,hessIncrement);		%+ B(h30,conj(q0))
  h31 = (2i*omega*eye(nphase)-jac)\(h31 - 6*g21*h20);  
  h22 = multilinear4(odefile,ten4,q0,q0,conj(q0),conj(q0),x,p,ten4Increment);		%  D(q0,q0,conj(q0),conj(q0))
  h22 = h22 + 4*multilinear3(odefile,tens,q0,conj(q0),h11,x,p,ten3Increment);		%+4C(q0,conj(q0),h11)
  h22 = h22 + 2*real(multilinear3(odefile,tens,conj(q0),conj(q0),h20,x,p,ten3Increment));	%+2*Re(C(q0,q0,h02))
  h22 = h22 + 4*real(multilinear2(odefile,hess,h21,conj(q0),x,p,hessIncrement));		%+2*Re(B(h21,conj(q0)))
  h22 = h22 + 2*multilinear2(odefile,hess,h11,h11,x,p,hessIncrement);			%+2B(h11,h11)
  h22 = -jac\(h22 + multilinear2(odefile,hess,h20,conj(h20),x,p,hessIncrement));	%+ B(h20,h02)
%5th order rhs
  h32 = multilinear5(odefile,ten5,q0,q0,q0,conj(q0),conj(q0),x,p,ten5Increment);	%  E(q0,q0,q0,conj(q0),conj(q0))
  h32 = h32 + 6*multilinear4(odefile,ten4,q0,q0,conj(q0),h11,x,p,ten4Increment);	%+6D(q0,q0,conj(q0),h11)
  h32 = h32 + 3*multilinear4(odefile,ten4,conj(q0),conj(q0),q0,h20,x,p,ten4Increment);	%+3D(conj(q0),conj(q0),q0,h20)
  h32 = h32 +   multilinear4(odefile,ten4,q0,q0,q0,conj(h20),x,p,ten4Increment);	%+ D(q0,q0,q0,h02)
  h32 = h32 + 6*multilinear3(odefile,tens,h11,h11,q0,x,p,ten3Increment);		%+6C(h11,h11,q0)
  h32 = h32 + 6*multilinear3(odefile,tens,conj(q0),h20,h11,x,p,ten3Increment);		%+6C(conj(q0),h20,h11)
  h32 = h32 + 6*multilinear3(odefile,tens,conj(q0),q0,h21,x,p,ten3Increment);		%+6C(conj(q0),q0,h21)
  h32 = h32 + 3*multilinear3(odefile,tens,q0,h20,conj(h20),x,p,ten3Increment);		%+3C(conj(q0),h20,h02)
  h32 = h32 + 3*multilinear3(odefile,tens,q0,q0,conj(h21),x,p,ten3Increment);		%+3C(q0,q0,h12)
  h32 = h32 +   multilinear3(odefile,tens,conj(q0),conj(q0),h30,x,p,ten3Increment);	%+ C(conj(q0),conj(q0),h30)
  h32 = h32 + 6*multilinear2(odefile,hess,h21,h11,x,p,hessIncrement);			%+6B(h21,h11)
  h32 = h32 + 3*multilinear2(odefile,hess,h22,q0,x,p,hessIncrement);			%+3B(h22,q0)
  h32 = h32 + 3*multilinear2(odefile,hess,h20,conj(h21),x,p,hessIncrement);		%+3B(h12,h20)
  h32 = h32 + 2*multilinear2(odefile,hess,h31,conj(q0),x,p,hessIncrement);		%+2B(h31,conj(q0))
  h32 = h32 +   multilinear2(odefile,hess,h30,conj(h20),x,p,hessIncrement);		%+ B(h30,h02)
  coef = real(p0'*h32/12);
