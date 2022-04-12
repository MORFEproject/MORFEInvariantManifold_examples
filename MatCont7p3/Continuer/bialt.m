function [D,N,R] = bialt(A)
%
% Computes indices of the bialternate product of 2A and I
%
% Syntax: [D,N,R] = bialt(A)

[n,n2] = size(A);
if n ~= n2
  msg=sprintf('bialt2ai: A must be square (was %d x %d)', n, n2);
  error(msg);
end

rn = n*(n-1)/2;
N = zeros (rn,rn);
R = zeros(rn+1,rn+1);
D = zeros(rn+1,rn+1);
j = 1:n-1; j = j.*(j-1)/2;

ps = 0;
for p = 2 : n
  qr = 1:p-1;

  R( ps+qr, p+j(p:(n-1)) ) = -A(qr,(p+1):n);

  qs = 0;
  for q = qr
    R( ps+q, [qs+(1:q-1) q+j(q:(n-1))] ) = [-A(p,1:q-1)  A(p,(q+1):n)];

    qs = qs+q-1;
  end %q
  R( ps+qr, ps+qr ) = A( qr, 1:p-1 );
  D(ps+qr,ps+qr) = eye(p-1)*A(p,p);
  ps = ps+p-1;
end %p
R = R(1:rn,1:rn);
D = D(1:rn,1:rn);
N = R; N(find(R>0)) = 0; N = -N;
R(find(R<0))=0;
