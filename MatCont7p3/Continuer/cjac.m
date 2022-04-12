function j=cjac(odefile,jacobian,x,p,ap)
%
% j=cjac(odefile,jacobian,x,p)
% Calculates the jacobian of the system. 
% p can be empty, if the active parameter is included in x (and the correct
% Jacobian function is passed).
% 
global cds

if nargin > 3 && ~isempty(p)
  
    fullx = [x;cell2mat(p(ap))];
    if (isfield(cds,'oldJacX')) && (length(cds.oldJacX) == length(fullx))
%     if isfield(cds,'oldJacX') && ~isempty(cds.oldJacX)
        % This checks if the Jacobian in this point has already been computed. If
        % so, the old one is recovered instead of recomputed.
        if fullx == cds.oldJacX & ~isempty(cds.oldJac)
            j = cds.oldJac;
            return;
        end
    end
    
    if (cds.options.SymDerivative >=1) 
        % Use symbolic derivatives if they are defined
        j = feval(jacobian, 0, x, p{:});
    else
        % If not, use finite differences
        nphase=size(x,1);
        for i=1: nphase
            x1 = x; x1(i) = x1(i)-cds.options.Increment;
            x2 = x; x2(i) = x2(i)+cds.options.Increment;
            j(:,i) = feval(odefile, 0, x2, p{:})-feval(odefile, 0, x1,  p{:});
        end
        j = j/(2*cds.options.Increment);
    end
    % Store the Jacobian globally for later use
    cds.oldJac = j;
    
else    
    fullx = x;
    if isfield(cds,'oldJacX') && (length(cds.oldJacX) == length(fullx)) 
%     if isfield(cds,'oldJacX') && ~isempty(cds.oldJacX)
        % This checks if the Jacobian in this point has already been computed. If
        % so, the old one is recovered instead of recomputed.
        if fullx == cds.oldJacX & ~isempty(cds.oldJacFull)
            j = cds.oldJacFull;
            return;
        end
    end
    if cds.symjac
        % Use curve-derivatives if they are defined
        j = feval(jacobian,x);
    else
        % If not, use finite differences
        x1 = x;
        x2 = x;
        for i=1:cds.ndim
            x1(i) = x1(i)-cds.options.Increment;
            x2(i) = x2(i)+cds.options.Increment;
            j(:,i) = feval(odefile, x2)-feval(odefile, x1);
            x1(i) = x(i);
            x2(i) = x(i);
        end
        j = j/(2*cds.options.Increment);
    end
    % Store the Jacobian globally for later use
    cds.oldJacFull = j;
end

% Store the Jacobian globally for later use
cds.oldJacX = fullx;