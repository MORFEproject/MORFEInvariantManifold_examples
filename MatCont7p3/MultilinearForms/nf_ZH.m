function coef = nf_ZH(odefile,jacobian,hessians,der3,x,p,nphase)
%
% coef = nf_ZH(odefile,jacobian,hessians,der3,x,p,nphase)
% compute normal form coefficients for zero-hopf.
%
global cds hds lpds
if ~isempty(hds)
  jac = cjac(odefile,jacobian,x,p,hds.ActiveParams);
else
  jac = cjac(odefile,jacobian,x,p,lpds.ActiveParams);
end    
  [X,D] = eig(jac);
  index1 = find( abs(diag(D))<1e-6 & sign(imag(diag(D)))==0);% fold part.
  index2 = find(abs(real(diag(D)))<1e-6 & sign(imag(diag(D)))==1);% hopf part.
if(isempty(index2)==1)
  debug('Neutral saddle\n');
  coef =[ 0,0,0];
  return;
end
  ev0 = diag(D(index1,index1));
  ev1 = diag(D(index2,index2));
  q0 = X(:,index1);
  q1 = X(:,index2);
  [X,D] = eig(jac');
  index1 = find( abs(diag(D))<1e-6 & sign(imag(diag(D)))== 0);% fold part.
  index2 = find(abs(real(diag(D)))<1e-6 & sign(imag(diag(D)))==-1);% hopf part.
  qad0 = X(:,index1);
  qad1 = X(:,index2);
  p0 = qad0/(q0'*qad0);
  p1 = qad1/(q1'*qad1);  
  hessIncrement = (cds.options.Increment)^(3.0/4.0);
  ten3Increment = (cds.options.Increment)^(3.0/5.0);
if (cds.options.SymDerivative >= 3)
  hess = chess(odefile,jacobian,hessians,x,p,ones(length(p),1));
  tens = ctens3(odefile,jacobian,hessians,der3,x,p,ones(length(p),1));
else
  hess = [];
  tens = [];
end
  h200 = multilinear2(odefile,hess,q0,q0,x,p,hessIncrement);		% B(q0,q0)
  h110 = multilinear2(odefile,hess,q0,q1,x,p,hessIncrement);		% B(q0,q1)
  h020 = multilinear2(odefile,hess,q1,q1,x,p,hessIncrement);		% B(q1,q1)
  h011 = multilinear2(odefile,hess,q1,conj(q1),x,p,hessIncrement);	% B(q1,conj(q1))
  f200 = p0'*h200/2.0; f011 = p0'*h011; g110 = p1'*h110;
  h200 = [jac q0 ; p0' 0] \ [2*f200*q0-h200 ; 0]; 			%-A^{INV}( B(q0,q0) - <p0,B(q0,q0)>q0 ) 
  h200 = h200(1:nphase);
  h110 = [(ev1*eye(nphase)-jac) q1 ; p1' 0] \ [h110-g110*q1 ; 0];
  h110 = h110(1:nphase);
  h020 = (2*ev1*eye(nphase) -jac)\h020;
  h011 = [jac q0 ; p0' 0] \ [f011*q0 - h011 ; 0];
  h011 = h011(1:nphase);
  h300 = multilinear3(odefile,tens,q0,q0,q0,x,p,ten3Increment);			%  C(q0,q0,q0)
  h300 = h300 +3*multilinear2(odefile,hess,q0,h200,x,p,hessIncrement);		%+3B(q0,h200)
  h111 = multilinear3(odefile,tens,q0,q1,conj(q1),x,p,ten3Increment);		%  C(q0,q1,conj(q1))
  h111 = h111 +2*real(multilinear2(odefile,hess,q1,conj(h110),x,p,hessIncrement));%+2Re(B(q1,h101))
  h111 = h111 + multilinear2(odefile,hess,q0,h011,x,p,hessIncrement);  		%+2B(q0,h011)
  h210 = multilinear3(odefile,tens,q0,q0,q1,x,p,ten3Increment); %  C(q0,q0,q1)
  h210 = h210 + 2*multilinear2(odefile,hess,q0,h110,x,p,hessIncrement); %+2B(q0,h110)
  h210 = h210 +   multilinear2(odefile,hess,q1,h200,x,p,hessIncrement); %+ B(q1,h200)
  h021 = multilinear3(odefile,tens,q1,q1,conj(q1),x,p,ten3Increment); %  C(q1,q1,conj(q1))
  h021 = h021 + 2*multilinear2(odefile,hess,q1,h011,x,p,hessIncrement); %+2B(q1,h011)
  h021 = h021 +   multilinear2(odefile,hess,conj(q1),h020,x,p,hessIncrement); %+ B(conj(q1),h020)
  f300 = p0'*h300/6.0; f111 = p0'*h111;
  g210 = p1'*h210/2.0; g021 = p1'*h021/2.0;
if(abs(f200*f011*real(g110)) < 1.0e-14)
  debug('Degenerate Zero-Hopf\n');
  coef =[ 0,0,0 ];
  return;
end
E0 = real(g210+g110*(real(g021)/f011-3*f300/f200/2+f111/2/f011)-g021*f200/f011);
coef = [ sign(f200*f011) real(g110/2/f200) sign(E0) ];
