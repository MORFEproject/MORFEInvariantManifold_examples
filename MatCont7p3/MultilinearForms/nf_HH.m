function coef = nf_HH(odefile,jacobian,hessians,der3,der4,der5,x,p,nphase)
%
% coef = nf_HH(odefile,jacobian,hessians,der3,x,p,nphase)
% compute normal form coefficients for double-hopf.
%
global cds hds
  jac=cjac(odefile,jacobian,x,p,hds.ActiveParams);
  [X,D] = eig(jac);
  D=diag(D);
  index=find(abs(real(D))<1e-6 & sign(imag(D))==1); % This should give a 1x2 vector
if(size(index)~=2)						% Otherwise a neutral saddle might be involved.
  debug('Neutral saddle ?\n');
  coef =[];
  return;
end
  if(imag(D(index(1))) < imag(D(index(2)))) 			% swap if necessary so that omega1>omega2
    index = [index(2);index(1)];
  end
  ev1 = D(index(1));
  ev2 = D(index(2));
  q1 = X(:,index(1));
  q2 = X(:,index(2));
  [XX,DD] = eig(jac');
  DD=diag(DD);
  index2=find(abs(real(DD))<1e-6 & sign(imag(DD))== -1); % This should again give a 1x2 vector
  if(imag(DD(index2(1))) > imag(DD(index2(2))))	 		% swap if necessary
    index2 = [index2(2);index2(1)];
  end
  qad1 = XX(:,index2(1));
  qad2 = XX(:,index2(2));
  p1=qad1/(q1'*qad1);
  p2=qad2/(q2'*qad2);
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
  ten4 = feval(hds.Der4,0,x, p{:});
  ten5 = feval(hds.Der5,0,x, p{:});
else
  ten4 = [];
  ten5 = [];
end
  Abor1 = [jac-ev1*eye(nphase) q1; p1' 0 ];
  Abor2 = [jac-ev2*eye(nphase) q2; p2' 0 ];
%2nd order vectors
  h2000 = (2*ev1*eye(nphase)-jac)\multilinear2(odefile,hess,q1,q1,x,p,hessIncrement);		% (2iw_1-A)\B(q1,q1)
  h1100 = -jac\multilinear2(odefile,hess,q1,conj(q1),x,p,hessIncrement);			% -A\B(q1,conj(q1))
  h1010 = ((ev1+ev2)*eye(nphase)-jac)\multilinear2(odefile,hess,q1,q2,x,p,hessIncrement);	% (i(w_1+w_2)-A)\B(q1,q2)
  h1001 = ((ev1-ev2)*eye(nphase)-jac)\multilinear2(odefile,hess,q1,conj(q2),x,p,hessIncrement);	% (i(w_1-w_2)-A)\B(q1,conj(q2))
  h0020 = (2*ev2*eye(nphase)-jac)\multilinear2(odefile,hess,q2,q2,x,p,hessIncrement);		% (2iw_2-A)\B(q2,q2)
  h0011 = -jac\multilinear2(odefile,hess,q2,conj(q2),x,p,hessIncrement);			% -A\B(q2,conj(q2))
  h0200 = conj(h2000); h0101 = conj(h1010); h0110 = conj(h1001); h0002 = conj(h0020);
%3rd order vectors
  h3000 = multilinear3(odefile,tens,q1,q1,q1,x,p,ten3Increment);				%  C(q1,q1,q1)
  h3000 = h3000 + 3*multilinear2(odefile,hess,q1,h2000,x,p,hessIncrement);			%+3B(h2000,q1)
  h3000 = (3*ev1*eye(nphase)-jac)\h3000; h0300 = conj(h3000);
  h2100 = multilinear3(odefile,tens,q1,q1,conj(q1),x,p,ten3Increment);				%  C(q1,q1,conj(q1))
  h2100 = h2100 + 2*multilinear2(odefile,hess,q1,h1100,x,p,hessIncrement);			%+2B(h1100,q1)
  h2100 = h2100 + multilinear2(odefile,hess,h2000,conj(q1),x,p,hessIncrement);			%+ B(h2000,conj(q1))  
  g2100 = p1'*h2100/2.0; h2100 = -Abor1\[2*g2100*q1-h2100;0];
  h2100 = h2100(1:nphase); h1200 = conj(h2100);
  h2010 = multilinear3(odefile,tens,q1,q1,q2,x,p,ten3Increment);				%  C(q1,q1,q2)
  h2010 = h2010 +2*multilinear2(odefile,hess,h1010,q1,x,p,hessIncrement);			%+2B(h1010,q1)
  h2010 = h2010 +  multilinear2(odefile,hess,h2000,q2,x,p,hessIncrement);			%+ B(h2000,q2)
  h2010 = ((2*ev1+ev2)*eye(nphase)-jac)\h2010; h0201 = conj(h2010);
  h2001 = multilinear3(odefile,tens,q1,q1,conj(q2),x,p,ten3Increment);				%  C(q1,q1,conj(q2))
  h2001 = h2001 +2*multilinear2(odefile,hess,h1001,q1,x,p,hessIncrement);			%+2B(h1001,q1)
  h2001 = h2001 +  multilinear2(odefile,hess,h2000,conj(q2),x,p,hessIncrement);			%+ B(h2000,conj(q2))
  h2001 = ((2*ev1-ev2)*eye(nphase)-jac)\h2001; h0210 = conj(h2001);
  h1110 = multilinear3(odefile,tens,q2,q1,conj(q1),x,p,ten3Increment);				%  C(q2,q1,conj(q1))
  h1110 = h1110 + multilinear2(odefile,hess,h1100,q2,x,p,hessIncrement);			%+ B(q2,h1100)
  h1110 = h1110 + multilinear2(odefile,hess,conj(h1001),q1,x,p,hessIncrement);			%+ B(q1,conj(h1001))
  h1110 = h1110 + multilinear2(odefile,hess,h1010,conj(q1),x,p,hessIncrement);			%+ B(conj(q1),h1010)
  g1110 = p2'*h1110; h1110 = -Abor2\[g1110*q2-h1110;0];
  h1110 = h1110(1:nphase); h1101 = conj(h1110);
  h1020 = multilinear3(odefile,tens,q2,q2,q1,x,p,ten3Increment);				%  C(q2,q2,q1)
  h1020 = h1020 +2*multilinear2(odefile,hess,h1010,q2,x,p,hessIncrement);			%+2B(h1010,q2)
  h1020 = h1020 +  multilinear2(odefile,hess,h0020,q1,x,p,hessIncrement);			%+ B(h0020,q1)
  h1020 = ((ev1+2*ev2)*eye(nphase)-jac)\h1020; h0102 = conj(h1020);
  h1011 = multilinear3(odefile,tens,q1,q2,conj(q2),x,p,ten3Increment);				%  C(q1,q2,conj(q2))
  h1011 = h1011 + multilinear2(odefile,hess,h0011,q1,x,p,hessIncrement);			%+ B(q1,h0011)
  h1011 = h1011 + multilinear2(odefile,hess,h1001,q2,x,p,hessIncrement);			%+ B(q2,h1001)
  h1011 = h1011 + multilinear2(odefile,hess,h1010,conj(q2),x,p,hessIncrement);			%+ B(conj(q2),h1010)
  g1011 = p1'*h1011; h1011 = -Abor1\[g1011*q1-h1011;0];
  h1011 = h1011(1:nphase); h0111 = conj(h1011);
  h1002 = multilinear3(odefile,tens,conj(q2),conj(q2),q1,x,p,ten3Increment);			%  C(conj(q2),conj(q2),q1)
  h1002 = h1002 +2*multilinear2(odefile,hess,h1001,conj(q2),x,p,hessIncrement);			%+2B(h1001,conj(q2))
  h1002 = h1002 +  multilinear2(odefile,hess,h0002,q1,x,p,hessIncrement);			%+ B(h0002,q1)
  h1002 = ((ev1-2*ev2)*eye(nphase)-jac)\h1002; h0120 = conj(h1002);
  h0030 = multilinear3(odefile,tens,q2,q2,q2,x,p,ten3Increment);				%  C(q2,q2,q2)
  h0030 = h0030 + 2*multilinear2(odefile,hess,h0020,q2,x,p,hessIncrement);			%+3B(h0020,q2)
  h0030 = (3*ev2*eye(nphase)-jac)\h0030; h0003 = conj(h0030);
  h0021 = multilinear3(odefile,tens,q2,q2,conj(q2),x,p,ten3Increment);				%  C(q2,q2,conj(q2))
  h0021 = h0021 + 2*multilinear2(odefile,hess,q2,h0011,x,p,hessIncrement);			%+2B(h0011,q2)
  h0021 = h0021 + multilinear2(odefile,hess,h0020,conj(q2),x,p,hessIncrement);			%+ B(h0020,conj(q2))
  g0021 = p2'*h0021/2.0; h0021 = -Abor2\[2*g0021*q2-h0021;0];
  h0021 = h0021(1:nphase); h0012 = conj(h0021);
%4rd order vectors
  h3100 = multilinear4(odefile,ten4,q1,q1,q1,conj(q1),x,p,ten4Increment);			%  D(q1,q1,q1,conj(q1))
  h3100 = h3100 + 3*multilinear3(odefile,tens,q1,q1,h1100,x,p,ten3Increment);			%+3C(q1,q1,h1100)
  h3100 = h3100 + 3*multilinear3(odefile,tens,q1,conj(q1),h2000,x,p,ten3Increment);		%+3C(q1,conj(q1),h2000)
  h3100 = h3100 + 3*multilinear2(odefile,hess,h3000,conj(q1),x,p,hessIncrement);		%+3B(h3000,conj(q1))
  h3100 = h3100 + 3*multilinear2(odefile,hess,h2100,q1,x,p,hessIncrement);			%+3B(h2100,q1)
  h3100 = (2*ev1*eye(nphase)-jac)\(h3100 - 6*g2100*h2000);
  h2200 = multilinear4(odefile,ten4,q1,q1,conj(q1),conj(q1),x,p,ten4Increment);			%  D(q1,q1,conj(q1),conj(q1))
  h2200 = h2200 + 4*multilinear3(odefile,tens,q1,conj(q1),h1100,x,p,ten3Increment);		%+4C(q1,conj(q1),h1100)
  h2200 = h2200 + 2*real(multilinear3(odefile,tens,q1,q1,h0200,x,p,ten3Increment));		%+ C(q1,q1,h0200)+ C(conj(q1),conj(q1),h2000)
  h2200 = h2200 + 4*real(multilinear2(odefile,hess,h2100,conj(q1),x,p,hessIncrement));		%+2B(h2100,conj(q1))+2B(h1200,q1)
  h2200 = h2200 + 2*multilinear2(odefile,hess,h1100,h1100,x,p,hessIncrement);			%+2B(h1100,h1100)
  h2200 = h2200 +   multilinear2(odefile,hess,h2000,h0200,x,p,hessIncrement);			%+ B(h2000,h0200)
  h2200 = -jac\(h2200 - 8*real(g2100)*h1100);
  h2110 = multilinear4(odefile,ten4,q1,q1,conj(q1),q2,x,p,ten4Increment);			%  D(q1,q1,conj(q1),q2)
  h2110 = h2110 + 2*multilinear3(odefile,tens,q1,conj(q1),h1010,x,p,ten3Increment);		%+2C(q1,conj(q1),h1010)
  h2110 = h2110 + 2*multilinear3(odefile,tens,q1,q2,h1100,x,p,ten3Increment);			%+2C(q1,q2,h1100)
  h2110 = h2110 +   multilinear3(odefile,tens,q1,q1,h0110,x,p,ten3Increment);			%+ C(q1,q1,h0110)
  h2110 = h2110 +   multilinear3(odefile,tens,q2,conj(q1),h2000,x,p,ten3Increment);		%+ C(q2,conj(q1),h2000)
  h2110 = h2110 + 2*multilinear2(odefile,hess,h1110,q1,x,p,hessIncrement);			%+2B(h1110,q1)
  h2110 = h2110 + 2*multilinear2(odefile,hess,h1100,h1010,x,p,hessIncrement);			%+2B(h1100,h1010)
  h2110 = h2110 +   multilinear2(odefile,hess,h2100,q2,x,p,hessIncrement);			%+ B(h2100,q2)
  h2110 = h2110 +   multilinear2(odefile,hess,h2000,h0110,x,p,hessIncrement);			%+ B(h2000,h0110)
  h2110 = h2110 +   multilinear2(odefile,hess,h2010,conj(q1),x,p,hessIncrement);		%+ B(h2010,conj(q1))
  h2110 = ((ev1+ev2)*eye(nphase)-jac)\(h2110 - 2*(g2100 + g1110)*h1010);
  h2101 = multilinear4(odefile,ten4,q1,q1,conj(q1),conj(q2),x,p,ten4Increment);			%  D(q1,q1,conj(q1),conj(q2))
  h2101 = h2101 + 2*multilinear3(odefile,tens,q1,conj(q1),h1001,x,p,ten3Increment);		%+2C(q1,conj(q1),h1001)
  h2101 = h2101 + 2*multilinear3(odefile,tens,q1,conj(q2),h1100,x,p,ten3Increment);		%+2C(q1,conj(q2),h1100)
  h2101 = h2101 +   multilinear3(odefile,tens,q1,q1,h0101,x,p,ten3Increment);			%+ C(q1,q1,h0101)
  h2101 = h2101 +   multilinear3(odefile,tens,conj(q2),conj(q1),h2000,x,p,ten3Increment);	%+ C(conj(q2),conj(q1),h2000)
  h2101 = h2101 + 2*multilinear2(odefile,hess,h1101,q1,x,p,hessIncrement);			%+2B(h1101,q1)
  h2101 = h2101 + 2*multilinear2(odefile,hess,h1100,h1001,x,p,hessIncrement);			%+2B(h1100,h1001)
  h2101 = h2101 +   multilinear2(odefile,hess,h2100,conj(q2),x,p,hessIncrement);		%+ B(h2100,conj(q2))
  h2101 = h2101 +   multilinear2(odefile,hess,h2000,h0101,x,p,hessIncrement);			%+ B(h2000,h0101)
  h2101 = h2101 +   multilinear2(odefile,hess,h2001,conj(q1),x,p,hessIncrement);		%+ B(h2001,conj(q1))
  h2101 = ((ev1-ev2)*eye(nphase)-jac)\(h2101 - 2*(g2100+conj(g1110))*h1001);h1210=conj(h2101);
  h2011 = multilinear4(odefile,ten4,q1,q1,q2,conj(q2),x,p,ten4Increment);			%  D(q1,q1,q2,conj(q2))
  h2011 = h2011 + 2*multilinear3(odefile,tens,q1,q2,h1001,x,p,ten3Increment);			%+2C(q1,q2,h1001)
  h2011 = h2011 + 2*multilinear3(odefile,tens,q1,conj(q2),h1010,x,p,ten3Increment);		%+2C(q1,conj(q2),h1010)
  h2011 = h2011 +   multilinear3(odefile,tens,q1,q1,h0011,x,p,ten3Increment);			%+ C(q1,q1,h0011)
  h2011 = h2011 +   multilinear3(odefile,tens,q2,conj(q2),h2000,x,p,ten3Increment);		%+ C(q2,conj(q2),h2000)
  h2011 = h2011 + 2*multilinear2(odefile,hess,h1010,h1001,x,p,hessIncrement);			%+2B(h1010,h1001)
  h2011 = h2011 +   multilinear2(odefile,hess,h2010,conj(q2),x,p,hessIncrement);		%+ B(h2010,conj(q2))
  h2011 = h2011 +   multilinear2(odefile,hess,h2001,q2,x,p,hessIncrement);			%+ B(h2001,q2)
  h2011 = h2011 +   multilinear2(odefile,hess,h2000,h0011,x,p,hessIncrement);			%+ B(h2000,h0011)
  h2011 = (2*ev1*eye(nphase)-jac)\(h2011 - 2*g1011*h2000);
  h1120 = multilinear4(odefile,ten4,q2,q2,q1,conj(q1),x,p,ten4Increment);			%  D(q2,q2,q1,conj(q1))
  h1120 = h1120 + 2*multilinear3(odefile,tens,q2,q1,h0110,x,p,ten3Increment);			%+2C(q2,q1,h0110)
  h1120 = h1120 + 2*multilinear3(odefile,tens,q2,conj(q1),h1010,x,p,ten3Increment);		%+2C(q2,conj(q1),h1010)
  h1120 = h1120 +   multilinear3(odefile,tens,q2,q2,h1100,x,p,ten3Increment);			%+ C(q2,q2,h1100)
  h1120 = h1120 +   multilinear3(odefile,tens,q1,conj(q1),h0020,x,p,ten3Increment);		%+ C(q1,conj(q1),h0020)
  h1120 = h1120 + 2*multilinear2(odefile,hess,h1010,h0110,x,p,hessIncrement);			%+2B(h1010,h0110)
  h1120 = h1120 +   multilinear2(odefile,hess,h1020,conj(q1),x,p,hessIncrement);		%+ B(h1020,conj(q1))
  h1120 = h1120 +   multilinear2(odefile,hess,h0120,q1,x,p,hessIncrement);			%+ B(h0120,q1)
  h1120 = h1120 +   multilinear2(odefile,hess,h0020,h1100,x,p,hessIncrement);			%+ B(h0020,h1100)
  h1120 = (2*ev2*eye(nphase)-jac)\(h1120-2*g1110*h0020);
  h1111 = multilinear4(odefile,ten4,q1,conj(q1),q2,conj(q2),x,p,ten4Increment);			%  D(q1,conj(q1),q2,conj(q2))
  h1111 = h1111 + 2*real(multilinear3(odefile,tens,q1,q2,h0101,x,p,ten3Increment));		%+2real(C(q1,q2,h0101))
  h1111 = h1111 + 2*real(multilinear3(odefile,tens,q1,conj(q2),h0110,x,p,ten3Increment));	%+2real(C(q1,conj(q2),h0110))
  h1111 = h1111 +   multilinear3(odefile,tens,q1,conj(q1),h0011,x,p,ten3Increment);		%+ C(q1,conj(q1),h0011)
  h1111 = h1111 +   multilinear3(odefile,tens,q2,conj(q2),h1100,x,p,ten3Increment);		%+ C(q2,conj(q2),h1100)
  h1111 = h1111 + 2*real(multilinear2(odefile,hess,q1,h0111,x,p,hessIncrement));		%+2real(B(q1,h0111))
  h1111 = h1111 + 2*real(multilinear2(odefile,hess,q2,h1101,x,p,hessIncrement));		%+2real(B(q2,h1101))
  h1111 = h1111 +   multilinear2(odefile,hess,h1001,h0110,x,p,hessIncrement);			%+ B(h1001,00110)
  h1111 = h1111 +   multilinear2(odefile,hess,h1010,h0101,x,p,hessIncrement);			%+ B(h1010,h0101)
  h1111 = -jac\(h1111-2*real(g1011)*h1100-2*real(g1110)*h0011);
  h1021 = multilinear4(odefile,ten4,q2,q2,conj(q2),q1,x,p,ten4Increment);			%  D(q2,q2,conj(q2),q1)
  h1021 = h1021 + 2*multilinear3(odefile,tens,q2,conj(q2),h1010,x,p,ten3Increment);		%+2C(q2,conj(q2),h1010)
  h1021 = h1021 + 2*multilinear3(odefile,tens,q2,q1,h0011,x,p,ten3Increment);			%+2C(q2,q1,h0011)
  h1021 = h1021 +   multilinear3(odefile,tens,q2,q2,h1001,x,p,ten3Increment);			%+ C(q2,q2,h1001)
  h1021 = h1021 +   multilinear3(odefile,tens,q1,conj(q2),h0020,x,p,ten3Increment);		%+ C(q1,conj(q2),h0020)
  h1021 = h1021 + 2*multilinear2(odefile,hess,h1011,q2,x,p,hessIncrement);			%+2B(h1011,q2)
  h1021 = h1021 + 2*multilinear2(odefile,hess,h0011,h1010,x,p,hessIncrement);			%+2B(h0011,h1010)
  h1021 = h1021 +   multilinear2(odefile,hess,h0021,q1,x,p,hessIncrement);			%+ B(h0021,q1)
  h1021 = h1021 +   multilinear2(odefile,hess,h0020,h1001,x,p,hessIncrement);			%+ B(h0020,h1001)
  h1021 = h1021 +   multilinear2(odefile,hess,h1020,conj(q2),x,p,hessIncrement);		%+ B(h1020,conj(q2))
  h1021 = ((ev1+ev2)*eye(nphase)-jac)\(h1021 - 2*(g0021 + g1011)*h1010);
  h1012 = multilinear4(odefile,ten4,q2,conj(q2),conj(q2),q1,x,p,ten4Increment);			%  D(q2,conj(q2),conj(q2),q1)
  h1012 = h1012 + 2*multilinear3(odefile,tens,q2,conj(q2),h1001,x,p,ten3Increment);		%+2C(q2,conj(q2),h1001)
  h1012 = h1012 + 2*multilinear3(odefile,tens,conj(q2),q1,h0011,x,p,ten3Increment);		%+2C(conj(q2),q1,h0011)
  h1012 = h1012 +   multilinear3(odefile,tens,conj(q2),conj(q2),h1010,x,p,ten3Increment);	%+ C(conj(q2),conj(q2),h1010)
  h1012 = h1012 +   multilinear3(odefile,tens,q1,q2,h0002,x,p,ten3Increment);			%+ C(q1,q2,h0002)
  h1012 = h1012 + 2*multilinear2(odefile,hess,h1011,conj(q2),x,p,hessIncrement);		%+2B(h1011,conj(q2))
  h1012 = h1012 + 2*multilinear2(odefile,hess,h0011,h1001,x,p,hessIncrement);			%+2B(h0011,h1001)
  h1012 = h1012 +   multilinear2(odefile,hess,h0012,q1,x,p,hessIncrement);			%+ B(h0012,q1)
  h1012 = h1012 +   multilinear2(odefile,hess,h0002,h1010,x,p,hessIncrement);			%+ B(h0002,h1010)
  h1012 = h1012 +   multilinear2(odefile,hess,h1002,q2,x,p,hessIncrement);			%+ B(h1002,q2)
  h1012 = ((ev1-ev2)*eye(nphase)-jac)\(h1012 - 2*(conj(g0021) + g1011)*h1001);h0121=conj(h1012);
  h0031 = multilinear4(odefile,ten4,q2,q2,q2,conj(q2),x,p,ten4Increment);			%  D(q2,q2,q2,conj(q2))
  h0031 = h0031 + 3*multilinear3(odefile,tens,q2,q2,h0011,x,p,ten3Increment);			%+3C(q2,q2,h0011)
  h0031 = h0031 + 3*multilinear3(odefile,tens,q2,conj(q2),h0020,x,p,ten3Increment);		%+3C(q2,conj(q2),h0020)
  h0031 = h0031 + 3*multilinear2(odefile,hess,h0030,conj(q2),x,p,hessIncrement);		%+3B(h0030,conj(q2))
  h0031 = h0031 + 3*multilinear2(odefile,hess,h0021,q2,x,p,hessIncrement);			%+3B(h0021,q2)
  h0031 = (2*ev1*eye(nphase)-jac)\(h0031 - 6*g0021*h0020);
  h0022 = multilinear4(odefile,ten4,q2,q2,conj(q2),conj(q2),x,p,ten4Increment);			%  D(q2,q2,conj(q2),conj(q2))
  h0022 = h0022 + 4*multilinear3(odefile,tens,q2,conj(q2),h0011,x,p,ten3Increment);		%+4C(q2,conj(q2),h0011)
  h0022 = h0022 + 2*real(multilinear3(odefile,tens,q2,q2,h0002,x,p,ten3Increment));		%+ C(q2,q2,h0002)+ C(conj(q2),conj(q2),h0020)
  h0022 = h0022 + 4*real(multilinear2(odefile,hess,h0021,conj(q2),x,p,hessIncrement));		%+2B(h2100,conj(q2))+2B(h0012,q2)
  h0022 = h0022 + 2*multilinear2(odefile,hess,h0011,h0011,x,p,hessIncrement);			%+2B(h0011,h0011)
  h0022 = h0022 +   multilinear2(odefile,hess,h0020,h0002,x,p,hessIncrement);			%+ B(h0020,h0002)
  h0022 = -jac\(h0022 - 8*real(g0021)*h0011);
%fifth order coefficients
%%%%%%%%%%%
  h3200 = multilinear5(odefile,ten5,q1,q1,q1,conj(q1),conj(q1),x,p,ten5Increment);		% E(q1,q1,q1,conj(q1),conj(q1))
  h3200 = h3200 + 6*multilinear4(odefile,ten4,q1,q1,conj(q1),h1100,x,p,ten4Increment);		%+6D(q1,q1,conj(q1),h1100)
  h3200 = h3200 + 3*multilinear4(odefile,ten4,conj(q1),conj(q1),q1,h2000,x,p,ten4Increment);	%+3D(conj(q1),conj(q1),q1,h2000)
  h3200 = h3200 +   multilinear4(odefile,ten4,q1,q1,q1,h0200,x,p,ten4Increment);		%+ D(q1,q1,q1,h0200)
  h3200 = h3200 + 6*multilinear3(odefile,tens,h1100,h1100,q1,x,p,ten3Increment);		%+6C(h1100,h1100,q1)
  h3200 = h3200 + 6*multilinear3(odefile,tens,conj(q1),h2000,h1100,x,p,ten3Increment);		%+6C(conj(q1),h2000,h1100)
  h3200 = h3200 + 6*multilinear3(odefile,tens,conj(q1),q1,h2100,x,p,ten3Increment);		%+6C(conj(q1),q1,h2100)
  h3200 = h3200 + 3*multilinear3(odefile,tens,q1,h2000,h0200,x,p,ten3Increment);		%+3C(q1,h2000,h0200)
  h3200 = h3200 + 3*multilinear3(odefile,tens,q1,q1,h1200,x,p,ten3Increment);			%+3C(q1,q1,h1200)
  h3200 = h3200 +   multilinear3(odefile,tens,conj(q1),conj(q1),h3000,x,p,ten3Increment);	%+ C(conj(q1),conj(q1),h3000)
  h3200 = h3200 + 6*multilinear2(odefile,hess,h2100,h1100,x,p,hessIncrement);			%+6B(h2100,h1100)
  h3200 = h3200 + 3*multilinear2(odefile,hess,h2200,q1,x,p,hessIncrement);			%+3B(h2200,q1)
  h3200 = h3200 + 3*multilinear2(odefile,hess,h2000,h1200,x,p,hessIncrement);			%+3B(h1200,h2000)
  h3200 = h3200 + 2*multilinear2(odefile,hess,h3100,conj(q1),x,p,hessIncrement);		%+2B(h3100,conj(q1))
  h3200 = h3200 +   multilinear2(odefile,hess,h3000,h0200,x,p,hessIncrement);			%+ B(h3000,h0200)
%%%%%%%%%%%
  h2111 = multilinear5(odefile,ten5,q1,q1,q2,conj(q1),conj(q2),x,p,ten5Increment);		% E(q1,q1,q2,conj(q1),conj(q2))
  h2111 = h2111 + 2*multilinear4(odefile,ten4,h1100,q1,q2,conj(q2),x,p,ten4Increment);		%+2D(h1100,q1,q2,conj(q2))
  h2111 = h2111 + 2*multilinear4(odefile,ten4,h1010,q1,conj(q1),conj(q2),x,p,ten4Increment);	%+2D(h1010,q1,conj(q1),conj(q2))
  h2111 = h2111 + 2*multilinear4(odefile,ten4,h1001,q1,q2,conj(q1),x,p,ten4Increment);		%+2D(h1001,q1,q2,conj(q1))
  h2111 = h2111 +   multilinear4(odefile,ten4,h0101,q1,q1,q2,x,p,ten4Increment);		%+ D(h0101,q1,q1,q2)
  h2111 = h2111 +   multilinear4(odefile,ten4,h0110,q1,q1,conj(q2),x,p,ten4Increment);		%+ D(h0110,q1,q1,conj(q2))
  h2111 = h2111 +   multilinear4(odefile,ten4,h2000,q2,conj(q1),conj(q2),x,p,ten4Increment);	%+ D(h2000,q2,conj(q1),conj(q2))
  h2111 = h2111 +   multilinear4(odefile,ten4,h0011,q1,q1,conj(q1),x,p,ten4Increment);		%+ D(h0011,q1,q1,conj(q1))
  h2111 = h2111 + 2*multilinear3(odefile,tens,h1101,q1,q2,x,p,ten3Increment);			%+2C(h1101,q1,q2)
  h2111 = h2111 + 2*multilinear3(odefile,tens,h1011,q1,conj(q1),x,p,ten3Increment);		%+2C(h1011,q1,conj(q1))
  h2111 = h2111 + 2*multilinear3(odefile,tens,h0110,h1001,q1,x,p,ten3Increment);		%+2C(h0110,h1001,q1)
  h2111 = h2111 + 2*multilinear3(odefile,tens,h1010,h1001,conj(q1),x,p,ten3Increment);		%+2C(h1010,h1001,conj(q1))
  h2111 = h2111 + 2*multilinear3(odefile,tens,h1010,h0101,q1,x,p,ten3Increment);		%+2C(h1010,h0101,q1)
  h2111 = h2111 + 2*multilinear3(odefile,tens,h1100,h0011,q1,x,p,ten3Increment);		%+2C(h1100,h0011,q1)
  h2111 = h2111 + 2*multilinear3(odefile,tens,h1100,h1001,q2,x,p,ten3Increment);		%+2C(h1100,h1001,q2)
  h2111 = h2111 + 2*multilinear3(odefile,tens,h1110,q1,conj(q2),x,p,ten3Increment);		%+2C(h1110,q1,conj(q2))
  h2111 = h2111 + 2*multilinear3(odefile,tens,h1100,h1010,conj(q2),x,p,ten3Increment);		%+2C(h1100,h1010,conj(q2))
  h2111 = h2111 +   multilinear3(odefile,tens,h2001,q2,conj(q1),x,p,ten3Increment);		%+ C(h2001,q2,conj(q1))
  h2111 = h2111 +   multilinear3(odefile,tens,h0111,q1,q1,x,p,ten3Increment);			%+ C(h0111,q1,q1)
  h2111 = h2111 +   multilinear3(odefile,tens,h2010,conj(q1),conj(q2),x,p,ten3Increment);	%+ C(h2010,conj(q1),conj(q2))
  h2111 = h2111 +   multilinear3(odefile,tens,h2100,q2,conj(q2),x,p,ten3Increment);		%+ C(h2100,q2,conj(q2))
  h2111 = h2111 +   multilinear3(odefile,tens,h2000,h0110,conj(q2),x,p,ten3Increment);		%+ C(h2000,h0110,conj(q2))
  h2111 = h2111 +   multilinear3(odefile,tens,h2000,h0011,conj(q1),x,p,ten3Increment);		%+ C(h2000,h0011,conj(q1))
  h2111 = h2111 +   multilinear3(odefile,tens,h2000,h0101,q2,x,p,ten3Increment);		%+ C(h2000,h0101,q2)
  h2111 = h2111 + 2*multilinear2(odefile,hess,h1010,h1101,x,p,hessIncrement);			%+2B(h1010,h1101)
  h2111 = h2111 + 2*multilinear2(odefile,hess,h1111,q1,x,p,hessIncrement);			%+2B(h1111,q1)
  h2111 = h2111 + 2*multilinear2(odefile,hess,h1110,h1001,x,p,hessIncrement);			%+2B(h1110,h1001)
  h2111 = h2111 + 2*multilinear2(odefile,hess,h1100,h1011,x,p,hessIncrement);			%+2B(h1100,h1011)
  h2111 = h2111 +   multilinear2(odefile,hess,h2011,conj(q1),x,p,hessIncrement);		%+ B(h2011,conj(q1))
  h2111 = h2111 +   multilinear2(odefile,hess,h2010,h0101,x,p,hessIncrement);			%+ B(h2010,h0101)
  h2111 = h2111 +   multilinear2(odefile,hess,h2101,q2,x,p,hessIncrement);			%+ B(h2101,q2)
  h2111 = h2111 +   multilinear2(odefile,hess,h2110,conj(q2),x,p,hessIncrement);		%+ B(h2110,conj(q2))
  h2111 = h2111 +   multilinear2(odefile,hess,h2100,h0011,x,p,hessIncrement);			%+ B(h2100,h0011)
  h2111 = h2111 +   multilinear2(odefile,hess,h0110,h2001,x,p,hessIncrement);			%+ B(h0110,h2001)
  h2111 = h2111 +   multilinear2(odefile,hess,h2000,h0111,x,p,hessIncrement);			%+ B(h2000,h0111)
%%%%%%%%%%%
  h1022 = multilinear5(odefile,ten5,q1,q2,q2,conj(q2),conj(q2),x,p,ten5Increment);		% E(q1,q2,q2,conj(q2),conj(q2))
  h1022 = h1022 + 4*multilinear4(odefile,ten4,h0011,q1,q2,conj(q2),x,p,ten4Increment);		%+4D(h0011,q1,q2,conj(q2))
  h1022 = h1022 + 2*multilinear4(odefile,ten4,h1001,q2,q2,conj(q2),x,p,ten4Increment);		%+2D(h1001,q2,q2,conj(q2))
  h1022 = h1022 + 2*multilinear4(odefile,ten4,h1010,q2,conj(q2),conj(q2),x,p,ten4Increment);	%+2D(h1010,q2,conj(q2),conj(q2))
  h1022 = h1022 +   multilinear4(odefile,ten4,h0002,q1,q2,q2,x,p,ten4Increment);		%+ D(h0002,q1,q2,q2)
  h1022 = h1022 +   multilinear4(odefile,ten4,h0020,q1,conj(q2),conj(q2),x,p,ten4Increment);	%+ D(h0020,q1,conj(q2),conj(q2))
  h1022 = h1022 + 4*multilinear3(odefile,tens,h1011,q2,conj(q2),x,p,ten3Increment);		%+4C(h1011,q2,conj(q2))
  h1022 = h1022 + 4*multilinear3(odefile,tens,h1001,h0011,q2,x,p,ten3Increment);		%+4C(h1001,h0011,q2)
  h1022 = h1022 + 4*multilinear3(odefile,tens,h1010,h0011,conj(q2),x,p,ten3Increment);		%+4C(h1010,h0011,conj(q2))
  h1022 = h1022 + 2*multilinear3(odefile,tens,h0012,q1,q2,x,p,ten3Increment);			%+2C(h0012,q1,q2)
  h1022 = h1022 + 2*multilinear3(odefile,tens,h0020,h1001,conj(q2),x,p,ten3Increment);		%+2C(h0020,h1001,conj(q2))
  h1022 = h1022 + 2*multilinear3(odefile,tens,h0021,q1,conj(q2),x,p,ten3Increment);		%+2C(h0021,q1,conj(q2))
  h1022 = h1022 + 2*multilinear3(odefile,tens,h1010,h0002,q2,x,p,ten3Increment);		%+2C(h1010,h0002,q2)
  h1022 = h1022 + 2*multilinear3(odefile,tens,h0011,h0011,q1,x,p,ten3Increment);		%+2C(h0011,h0011,q1)
  h1022 = h1022 +   multilinear3(odefile,tens,h1002,q2,q2,x,p,ten3Increment);			%+ C(h1002,q2,q2)
  h1022 = h1022 +   multilinear3(odefile,tens,h1020,conj(q2),conj(q2),x,p,ten3Increment);	%+ C(h1020,conj(q2),conj(q2))
  h1022 = h1022 +   multilinear3(odefile,tens,h0020,h0002,q1,x,p,ten3Increment);		%+ C(h0020,h0002,q1)
  h1022 = h1022 + 4*multilinear2(odefile,hess,h0011,h1011,x,p,hessIncrement);			%+4B(h0011,h1011)
  h1022 = h1022 + 2*multilinear2(odefile,hess,h1021,conj(q2),x,p,hessIncrement);		%+2B(h1021,conj(q2))
  h1022 = h1022 + 2*multilinear2(odefile,hess,h1001,h0021,x,p,hessIncrement);			%+2B(h1001,h0021)
  h1022 = h1022 + 2*multilinear2(odefile,hess,h1012,q2,x,p,hessIncrement);			%+2B(h1012,q2)
  h1022 = h1022 + 2*multilinear2(odefile,hess,h1010,h0012,x,p,hessIncrement);			%+2B(h1010,h0012)
  h1022 = h1022 +   multilinear2(odefile,hess,q1,h0022,x,p,hessIncrement);			%+ B(q1,h0022)
  h1022 = h1022 +   multilinear2(odefile,hess,h1020,h0002,x,p,hessIncrement);			%+ B(h1020,h0002)
  h1022 = h1022 +   multilinear2(odefile,hess,h0020,h1002,x,p,hessIncrement);			%+ B(h0020,h1002)
%%%%%%%%%%%
  h2210 = multilinear5(odefile,ten5,q2,q1,q1,conj(q1),conj(q1),x,p,ten5Increment);		% E(q2,q1,q1,conj(q1),conj(q1))
  h2210 = h2210 + 4*multilinear4(odefile,ten4,h1100,q1,conj(q1),q2,x,p,ten4Increment);		%+4D(h1100,q1,conj(q1),q2)
  h2210 = h2210 + 2*multilinear4(odefile,ten4,h1010,q1,conj(q1),conj(q1),x,p,ten4Increment);	%+2D(h1010,q1,conj(q1),conj(q1))
  h2210 = h2210 + 2*multilinear4(odefile,ten4,h0110,q1,q1,conj(q1),x,p,ten4Increment);		%+2D(h0110,q1,q1,conj(q1))
  h2210 = h2210 +   multilinear4(odefile,ten4,h2000,conj(q1),conj(q1),q2,x,p,ten4Increment);	%+ D(h2000,conj(q1),conj(q1),q2)
  h2210 = h2210 +   multilinear4(odefile,ten4,h0200,q1,q1,q2,x,p,ten4Increment);		%+ D(h0200,q1,q1,q2)
  h2210 = h2210 + 4*multilinear3(odefile,tens,h1110,q1,conj(q1),x,p,ten3Increment);		%+4C(h1110,q1,conj(q1))
  h2210 = h2210 + 4*multilinear3(odefile,tens,h1100,h0110,q1,x,p,ten3Increment);		%+4C(h1100,h0110,q1)
  h2210 = h2210 + 4*multilinear3(odefile,tens,h1100,h1010,conj(q1),x,p,ten3Increment);		%+4C(h1100,h1010,conj(q1))
  h2210 = h2210 + 2*multilinear3(odefile,tens,h1100,h1100,q2,x,p,ten3Increment);		%+2C(h1100,h1100,q2)
  h2210 = h2210 + 2*multilinear3(odefile,tens,h1200,q1,q2,x,p,ten3Increment);			%+2C(h1200,q1,q2)
  h2210 = h2210 + 2*multilinear3(odefile,tens,h2100,conj(q1),q2,x,p,ten3Increment);		%+2C(h2100,conj(q1),q2)
  h2210 = h2210 + 2*multilinear3(odefile,tens,h0200,h1010,q1,x,p,ten3Increment);		%+2C(h0200,h1010,q1)
  h2210 = h2210 + 2*multilinear3(odefile,tens,h2000,h0110,conj(q1),x,p,ten3Increment);		%+2C(h2000,h0110,conj(q1))
  h2210 = h2210 +   multilinear3(odefile,tens,h2000,h0200,q2,x,p,ten3Increment);		%+ C(h2000,h0200,q2)
  h2210 = h2210 +   multilinear3(odefile,tens,h0210,q1,q1,x,p,ten3Increment);			%+ C(h0210,q1,q1)
  h2210 = h2210 +   multilinear3(odefile,tens,h2010,conj(q1),conj(q1),x,p,ten3Increment);	%+ C(h2010,conj(q1),conj(q1))
  h2210 = h2210 + 4*multilinear2(odefile,hess,h1100,h1110,x,p,hessIncrement);			%+4B(h1100,h1110)
  h2210 = h2210 + 2*multilinear2(odefile,hess,h1210,q1,x,p,hessIncrement);			%+2B(h1210,q1)
  h2210 = h2210 + 2*multilinear2(odefile,hess,h2110,conj(q1),x,p,hessIncrement);		%+2B(h2110,conj(q1))
  h2210 = h2210 + 2*multilinear2(odefile,hess,h2100,h0110,x,p,hessIncrement);			%+2B(h2100,h0110)
  h2210 = h2210 + 2*multilinear2(odefile,hess,h1200,h1010,x,p,hessIncrement);			%+2B(h1200,h1010)
  h2210 = h2210 +   multilinear2(odefile,hess,h2200,q2,x,p,hessIncrement);			%+ B(h2200,q2)
  h2210 = h2210 +   multilinear2(odefile,hess,h2000,h0210,x,p,hessIncrement);			%+ B(h2000,h0210)
  h2210 = h2210 +   multilinear2(odefile,hess,h0200,h2010,x,p,hessIncrement);			%+ B(h0200,h2010)
%%%%%%%%%%%
  h1121 = multilinear5(odefile,ten5,q1,q2,q2,conj(q1),conj(q2),x,p,ten5Increment);		% E(q1,q2,q2,conj(q1),conj(q2))
  h1121 = h1121 + 2*multilinear4(odefile,ten4,h1010,conj(q1),q2,conj(q2),x,p,ten4Increment);	%+2D(h1010,conj(q1),q2,conj(q2))
  h1121 = h1121 + 2*multilinear4(odefile,ten4,h0110,q1,q2,conj(q2),x,p,ten4Increment);		%+2D(h0110,q1,q2,conj(q2))
  h1121 = h1121 + 2*multilinear4(odefile,ten4,h0011,q1,q2,conj(q1),x,p,ten4Increment);		%+2D(h0011,q1,q2,conj(q1))
  h1121 = h1121 +   multilinear4(odefile,ten4,h0101,q1,q2,q2,x,p,ten4Increment);		%+ D(h0101,q1,q2,q2)
  h1121 = h1121 +   multilinear4(odefile,ten4,h1001,conj(q1),q2,q2,x,p,ten4Increment);		%+ D(h1001,conj(q1),q2,q2)
  h1121 = h1121 +   multilinear4(odefile,ten4,h1100,q2,q2,conj(q2),x,p,ten4Increment);		%+ D(h1100,q2,q2,conj(q2))
  h1121 = h1121 +   multilinear4(odefile,ten4,h0020,q1,conj(q1),conj(q2),x,p,ten4Increment);	%+ D(h0020,q1,conj(q1),conj(q2))
  h1121 = h1121 + 2*multilinear3(odefile,tens,h1011,conj(q1),q2,x,p,ten3Increment);		%+2C(h1011,conj(q1),q2)
  h1121 = h1121 + 2*multilinear3(odefile,tens,h0111,q1,q2,x,p,ten3Increment);			%+2C(h0111,q1,q2)
  h1121 = h1121 + 2*multilinear3(odefile,tens,h1110,q2,conj(q2),x,p,ten3Increment);		%+2C(h1110,q2,conj(q2))
  h1121 = h1121 + 2*multilinear3(odefile,tens,h0110,h1001,q2,x,p,ten3Increment);		%+2C(h0110,h1001,q2)
  h1121 = h1121 + 2*multilinear3(odefile,tens,h1100,h0011,q2,x,p,ten3Increment);		%+2C(h1100,h0011,q2)
  h1121 = h1121 + 2*multilinear3(odefile,tens,h1010,h0101,q2,x,p,ten3Increment);		%+2C(h1010,h0101,q2)
  h1121 = h1121 + 2*multilinear3(odefile,tens,h1010,h0110,conj(q2),x,p,ten3Increment);		%+2C(h1010,h0110,conj(q2))
  h1121 = h1121 + 2*multilinear3(odefile,tens,h1010,h0011,conj(q1),x,p,ten3Increment);		%+2C(h1010,h0011,conj(q1))
  h1121 = h1121 + 2*multilinear3(odefile,tens,h0110,h0011,q1,x,p,ten3Increment);		%+2C(h0110,h0011,q1)
  h1121 = h1121 +   multilinear3(odefile,tens,h1020,conj(q1),conj(q2),x,p,ten3Increment);	%+ C(h1020,conj(q1),conj(q2))
  h1121 = h1121 +   multilinear3(odefile,tens,h0021,q1,conj(q1),x,p,ten3Increment);		%+ C(h0021,q1,conj(q1))
  h1121 = h1121 +   multilinear3(odefile,tens,h0120,q1,conj(q2),x,p,ten3Increment);		%+ C(h0120,q1,conj(q2))
  h1121 = h1121 +   multilinear3(odefile,tens,h1101,q2,q2,x,p,ten3Increment);			%+ C(h1101,q2,q2)
  h1121 = h1121 +   multilinear3(odefile,tens,h0020,h0101,q1,x,p,ten3Increment);		%+ C(h0020,h0101,q1)
  h1121 = h1121 +   multilinear3(odefile,tens,h0020,h1001,conj(q1),x,p,ten3Increment);		%+ C(h0020,h1001,conj(q1))
  h1121 = h1121 +   multilinear3(odefile,tens,h0110,h0011,q1,x,p,ten3Increment);		%+ C(h0110,h0011,q1)
  h1121 = h1121 + 2*multilinear2(odefile,hess,h1010,h0111,x,p,hessIncrement);			%+2B(h1010,h0111)
  h1121 = h1121 + 2*multilinear2(odefile,hess,h1111,q2,x,p,hessIncrement);			%+2B(h1111,q2)
  h1121 = h1121 + 2*multilinear2(odefile,hess,h0110,h1011,x,p,hessIncrement);			%+2B(h0110,h1011)
  h1121 = h1121 + 2*multilinear2(odefile,hess,h0110,h1011,x,p,hessIncrement);			%+2B(h0110,h1011)
  h1121 = h1121 + 2*multilinear2(odefile,hess,h1110,h0011,x,p,hessIncrement);			%+2B(h1110,h0011)
  h1121 = h1121 +   multilinear2(odefile,hess,h0121,q1,x,p,hessIncrement);			%+ B(h0121,q1)
  h1121 = h1121 +   multilinear2(odefile,hess,h1021,conj(q1),x,p,hessIncrement);		%+ B(h1021,conj(q1))
  h1121 = h1121 +   multilinear2(odefile,hess,h1120,conj(q2),x,p,hessIncrement);		%+ B(h1120,conj(q2))
  h1121 = h1121 +   multilinear2(odefile,hess,h1020,h0101,x,p,hessIncrement);			%+ B(h1020,h0101)
  h1121 = h1121 +   multilinear2(odefile,hess,h0020,h1101,x,p,hessIncrement);			%+ B(h0020,h1101)
  h1121 = h1121 +   multilinear2(odefile,hess,h1100,h0021,x,p,hessIncrement);			%+ B(h1100,h0021)
  h1121 = h1121 +   multilinear2(odefile,hess,h0120,h1001,x,p,hessIncrement);			%+ B(h0120,h1001)
%%%%%%%%%%%
  h0032 = multilinear5(odefile,ten5,q2,q2,q2,conj(q2),conj(q2),x,p,ten5Increment);		%  E(q2,q2,q2,conj(q2),conj(q2))
  h0032 = h0032 + 6*multilinear4(odefile,ten4,q2,q2,conj(q2),h0011,x,p,ten4Increment);		%+6D(q2,q2,conj(q2),h0011)
  h0032 = h0032 + 3*multilinear4(odefile,ten4,conj(q2),conj(q2),q2,h0020,x,p,ten4Increment);	%+3D(conj(q2),conj(q2),q2,h0020)
  h0032 = h0032 +   multilinear4(odefile,ten4,q2,q2,q2,h0002,x,p,ten4Increment);		%+ D(q2,q2,q2,h0002)
  h0032 = h0032 + 6*multilinear3(odefile,tens,h0011,h0011,q2,x,p,ten3Increment);		%+6C(h0011,h0011,q2)
  h0032 = h0032 + 6*multilinear3(odefile,tens,conj(q2),h0020,h0011,x,p,ten3Increment);		%+6C(conj(q2),h0020,h0011)
  h0032 = h0032 + 6*multilinear3(odefile,tens,conj(q2),q2,h0021,x,p,ten3Increment);		%+6C(conj(q2),q2,h0021)
  h0032 = h0032 + 3*multilinear3(odefile,tens,q2,h0020,h0002,x,p,ten3Increment);		%+3C(conj(q2),h0020,h0002)
  h0032 = h0032 + 3*multilinear3(odefile,tens,q2,q2,h0012,x,p,ten3Increment);			%+3C(q2,q2,h0012)
  h0032 = h0032 +   multilinear3(odefile,tens,conj(q2),conj(q2),h0030,x,p,ten3Increment);	%+ C(conj(q2),conj(q2),h0030)
  h0032 = h0032 + 6*multilinear2(odefile,hess,h0021,h0011,x,p,hessIncrement);			%+6B(h0021,h0011)
  h0032 = h0032 + 3*multilinear2(odefile,hess,h0022,q2,x,p,hessIncrement);			%+3B(h0022,q2)
  h0032 = h0032 + 3*multilinear2(odefile,hess,h0020,h0012,x,p,hessIncrement);			%+3B(h0020,h0012)
  h0032 = h0032 + 2*multilinear2(odefile,hess,h0031,conj(q2),x,p,hessIncrement);		%+2B(h0031,conj(q2))
  h0032 = h0032 +   multilinear2(odefile,hess,h0030,h0002,x,p,hessIncrement);			%+ B(h0030,h0002)
%------------%
%coefficients%
%------------%
  g2100 = real(g2100); g1011 = real(g1011); g1110 = real(g1110); g0021 = real(g0021);
  g3200 = real(p1'*h3200/12); g2111 = real(p1'*h2111/2); g1022 = real(p1'*h1022/4);
  g0032 = real(p2'*h0032/12); g1121 = real(p2'*h1121/2); g2210 = real(p2'*h2210/4);
  s1 = g1022 + g1011*(g1121/g1110 - 2*g0032/g0021 - g3200*g0021/g2100/g1110);
  s2 = g2210 + g1110*(g2111/g1011 - 2*g3200/g2100 - g0032*g2100/g0021/g1011);
  coef = [ sign(g2100*g0021) g1011/g0021 g1110/g2100 s1/g0021/g0021 s2/g2100/g2100];
