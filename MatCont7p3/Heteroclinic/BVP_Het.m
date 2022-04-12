function f = BVP_Het(x,x0,x1,p,T,eps0,eps1,YS,YU)

global hetds 

% extract ups
ups = reshape(x,hetds.nphase,hetds.tps);
p = num2cell(p);

% Cycle-component
% ---------------
% Here the cycle is approximated in each mesh interval as a polynomial
% of degree nCol, using the Lagrange polynomials. We approximate using the
% collocation points in the subfunction below, and here we compute the
% values of the polynomials at the mesh points.
range1 = hetds.cols_p1;
range2 = hetds.phases;
f = [];
for j=hetds.tsts
    % value of polynomial on each collocation point
    xval = ups(:,range1)*hetds.wt;  
    % derivative of polynomial on each collocation point
    derval  = ups(:,range1)*hetds.wpvec/hetds.dt(j);
    % evaluate function value on each collocation point
    for c=hetds.cols
        f(range2,1) = derval(:,c) - 2*T * feval(hetds.func, 0, xval(:,c), p{:});        
        % Shift for storage
        range2 = range2+hetds.nphase;
    end
    range1 = range1+hetds.ncol;
end

% Component 2
% -----------
% equilibrium condition
f(end+1:end+hetds.nphase,1) = feval(hetds.func, 0, x0, p{:}); 
f(end+1:end+hetds.nphase,1) = feval(hetds.func, 0, x1, p{:}); 
  
% Component 3
% -----------
if sum(hetds.extravec) == 2
    % Integral constraint (phase condition)
    if isempty(hetds.upold)
        hetds.upold = ups;
    end
    dotprod = dot(ups-hetds.upold, hetds.upoldp);
    dotprodmatrix = dotprod(hetds.idxmat);
    f(end+1,1) = sum(hetds.dt .* (hetds.wi * dotprodmatrix));
elseif sum(hetds.extravec) ~= 1
    error('Wrong number of free heteroclinic parameters');
end

% Component 4
% -----------
% Ricatti equations
if (hetds.nphase-hetds.npos)
    f(end+1:end+(hetds.nphase-hetds.npos)*hetds.npos) = RicattiEvalHet(x0,p,1,YU);
end
if (hetds.nphase-hetds.nneg)
    f(end+1:end+(hetds.nphase-hetds.nneg)*hetds.nneg) = RicattiEvalHet(x1,p,0,YS);    
end

% Component 5
% -----------
% Last and first vectors along stable and unstable eigenspaces
Q0S = hetds.oldStableQ;
Q0U = hetds.oldUnstableQ;
if (hetds.nphase-hetds.npos)
    Q1U = Q0U * [-YU'; eye(size(YU,1))];
    for i=1:hetds.nphase-hetds.npos
        f(end+1,1) = (ups(:,1) - x0)' * Q1U(:,end-i+1);
    end
end
if (hetds.nphase-hetds.nneg)
    Q1S = Q0S * [-YS'; eye(size(YS,1))];
    for i=1:hetds.nphase-hetds.nneg
        f(end+1,1) = (ups(:,end) - x1)' * Q1S(:,end-i+1);       
    end
end

% Component 5
% -----------
% Distances from endpoints to equilibrium equal to epsilons
f(end+1,1) = norm(ups(:,1) - x0) - eps0;
f(end+1,1) = norm(ups(:,end) - x1) - eps1;

