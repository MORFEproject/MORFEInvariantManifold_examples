function K = fastkron(c,p,A,B)

t = p:((c+2)*p-1);
K = A(ones(1,p),fix(t/p)).*B(:,rem(t,p)+1);
