function j=cjacp(odefile,jacobianp,x,p,ap)
%
% j=cjacp(odefile,jacobianp,x,p,ap)
% Calculates the jacobian of the system with respect to the active
% parameter ap.
% 

global cds
if cds.options.SymDerivativeP >= 1
    % Use symbolic derivatives if they are defined
    j = feval(jacobianp, 0, x, p{:});
else
    % If not, use finite differences
    if size(ap,1) ~= 1
        ap = ap';
    end
    for i=ap
        p1 = p; p1{i} = p1{i}-cds.options.Increment;
        p2 = p; p2{i} = p2{i}+cds.options.Increment;
        j(:,i) = feval(odefile, 0, x, p2{:})-feval(odefile, 0, x, p1{:});
    end
    j = j/(2*cds.options.Increment);
end
j=j(:,ap);
