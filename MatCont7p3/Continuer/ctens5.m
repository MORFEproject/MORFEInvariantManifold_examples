function h=ctens5(odefile,jacobian,hessians,tensor3,tensor4,tensor5,x,p,ap)
global  cds

nphase = size(x,1);
if cds.options.SymDerivative >= 5
    h = feval(tensor5, 0, x, p{:});
else
    for i=1:nphase
        x1 = x; x1(i) = x1(i)-cds.options.Increment;
        x2 = x; x2(i) = x2(i)+cds.options.Increment;
        h(:,:,:,:,:,i) = ctens4(odefile,jacobian,hessians,tensor3,tensor4,x2,p,ap)-ctens4(odefile,jacobian,hessians,tensor3,tensor4,x1,p,ap);
    end
    h = h/(2*cds.options.Increment);
end
