function h=chess(odefile,jacobian,hessians,x,p,ap)
%
% h=chess(odefile,jacobian,hessians,x,p)
% Calculates the hessian of the system
% 
global cds

if nargin > 4 && ~isempty(p)
    if cds.options.SymDerivative >= 2
        % Use symbolic derivatives if they are defined
        h = feval(hessians, 0, x, p{:});
    else
        % If not, use finite differences
        nphase=size(x,1);
        for i=1:nphase
            x1 = x; x1(i) = x1(i)-cds.options.Increment;
            x2 = x; x2(i) = x2(i)+cds.options.Increment;
            h(:,:,i) = cjac(odefile,jacobian,x2,p,ap)-cjac(odefile,jacobian,x1,p,ap);
        end
        h = h/(2*cds.options.Increment);
    end
    
else
    if cds.symhess
        % Use symbolic derivatives if they are defined
        h = feval(hessians, x);
    else
        % If not, use finite differences
        for i=1:cds.ndim
            x1 = x; x1(i) = x1(i)-cds.options.Increment;
            x2 = x; x2(i) = x2(i)+cds.options.Increment;
            h(:,:,i) = cjac(odefile,jacobian,x2,[])-cjac(odefile,jacobian,x1,[]);
        end
        h = h/(2*cds.options.Increment);
    end
end
    
