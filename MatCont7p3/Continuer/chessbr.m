function h=chessbr(odefile,jacobian,hessiansp,x,p,ap,bp)
%
% h=chessp(odefile,jacobianp,hessiansp,x,p,ap)
% Calculates the hessian of the system with respect to the active
% parameter ap.
% 
global cds

if cds.options.SymDerivativeP >= 2
    % Use symbolic derivatives if they are defined
    h = feval(hessiansp, 0, x, p{:});
    h = h(:,:,bp);
else
    % If not, use finite differences
    if size(bp,1) ~= 1
        bp = bp';
    end
    for i=bp
        p1 = p; p1{i} = p1{i}-cds.options.Increment;
        p2 = p; p2{i} = p2{i}+cds.options.Increment;
        h(:,:,i) = cjac(odefile,jacobian,x,p2,ap)-cjac(odefile,jacobian,x,p1,ap);
    end
    h = h(:,:,bp)/(2*cds.options.Increment);
end
