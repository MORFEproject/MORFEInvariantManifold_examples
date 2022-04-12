function h=chesspbr(odefile,jacobianp,hessiansp,x,p,ap,bp)
%
% h=chessp(odefile,jacobianp,hessiansp,x,p,ap)
% Calculates the hessian of the system with respect to the active
% parameter ap.
% 
global cds

    % If not, use finite differences
    if size(ap,1) ~= 1
        ap = ap';
    end
    for i=ap
        p1 = p; p1{i} = p1{i}-cds.options.Increment;
        p2 = p; p2{i} = p2{i}+cds.options.Increment;
        h(:,:,i) = cjacbr(odefile,jacobianp,x,p2,ap,bp)-cjacbr(odefile,jacobianp,x,p1,ap,bp);
    end
    h = h(:,:,ap)/(2*cds.options.Increment);
