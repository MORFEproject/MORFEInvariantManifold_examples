function h=ctens3(odefile,jacobian,hessians,tensor3,x,p,ap)
global  cds

nphase = size(x,1);
if cds.options.SymDerivative >= 3
    h = feval(tensor3, 0, x, p{:});
else
    for i=1:nphase
        x1 = x; x1(i) = x1(i)-cds.options.Increment;
        x2 = x; x2(i) = x2(i)+cds.options.Increment;
        h(:,:,:,i) = chess(odefile,jacobian,hessians,x2,p,ap)-chess(odefile,jacobian,hessians,x1,p,ap);
    end
    h = h/(2*cds.options.Increment);
end
