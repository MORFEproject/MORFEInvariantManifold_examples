function f = BVP_Hom(x,x0,p,T,eps0,eps1,YS,YU)
global homds 

% extract ups
ups = reshape(x,homds.nphase,homds.tps);
p = num2cell(p);

% Cycle-component
% ---------------
% Here the cycle is approximated in each mesh interval as a polynomial
% of degree nCol, using the Lagrange polynomials. We approximate using the
% collocation points in the subfunction below, and here we compute the
% values of the polynomials at the mesh points.
range1 = homds.cols_p1;
range2 = homds.phases;
f = [];
for j=homds.tsts
    % value of polynomial on each collocation point
    xval = ups(:,range1)*homds.wt;  
    % derivative of polynomial on each collocation point
    derval  = ups(:,range1)*homds.wpvec/homds.dt(j);
    % evaluate function value on each collocation point
    for c=homds.cols
        f(range2,1) = derval(:,c) - 2*T * feval(homds.func, 0, xval(:,c), p{:});        
        % Shift for storage
        range2 = range2+homds.nphase;
    end
    range1 = range1+homds.ncol;
end

% Component 2
% -----------
% equilibrium condition
f(end+1:end+homds.nphase,1) = feval(homds.func, 0, x0, p{:}); 

% Component 3
% -----------
if sum(homds.extravec) == 2
    % Integral constraint (phase condition)
    if isempty(homds.upold)
        homds.upold = ups;
    end
    dotprod = dot(ups-homds.upold, homds.upoldp);
    dotprodmatrix = dotprod(homds.idxmat);
    f(end+1,1) = sum(homds.dt .* (homds.wi * dotprodmatrix));
elseif sum(homds.extravec) ~= 1
    error('Wrong number of free homoclinic parameters');
end

% Component 4
% -----------
% Ricatti equations
f(end+1:end+2*homds.nneg*homds.npos) = ricattiEval(x0,p,YU,YS);

% Component 5
% -----------
% Last and first vectors along stable and unstable eigenspaces
Q0S = homds.oldStableQ;
Q0U = homds.oldUnstableQ;

if homds.npos
    Q1U = Q0U * [-YU'; eye(size(YU,1))];
    for i=1:homds.nneg
        f(end+1,1) = (ups(:,1) - x0)' * Q1U(:,end-i+1);
    end
end
if homds.nneg
    Q1S = Q0S * [-YS'; eye(size(YS,1))];
    for i=1:homds.npos
        f(end+1,1) = (ups(:,end) - x0)' * Q1S(:,end-i+1);
    end
end

% Component 6
% -----------
% Distances from endpoints to equilibrium equal to epsilons
f(end+1,1) = norm(ups(:,1) - x0) - eps0;
f(end+1,1) = norm(ups(:,end) - x0) - eps1;

