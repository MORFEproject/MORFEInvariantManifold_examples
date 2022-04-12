function f = BVP_HSN(x,x0,p,T,eps0,eps1,YS,YU)
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
    
% Component 2bis
% --------------
% limit point (saddle-node) condition
jac = cjac(homds.func,homds.Jacobian,x0,p,homds.ActiveParams);
Bord = [   jac       homds.wvector;...
       homds.vvector'      0      ];
bunit = [zeros(homds.nphase,1);1];
vext = Bord \ bunit;
f(end+1,1) = vext(end);    

% Component 3
% -----------
% Integral constraint (phase condition)
if isempty(homds.upold)
    homds.upold = ups;
end
dotprod = dot(ups-homds.upold, homds.upoldp);
dotprodmatrix = dotprod(homds.idxmat);
f(end+1,1) = sum(homds.dt .* (homds.wi * dotprodmatrix));

% Component 4
% -----------
% Ricatti equations
tmp = HSN_ricattiEval(x0,p,YU,YS);
f(end+1:end+length(tmp)) = tmp;

% Component 5
% -----------
% Last and first vectors along stable and unstable eigenspaces
Q0S = homds.oldStableQ;
Q0U = homds.oldUnstableQ;
opts.disp = 0;
[V,D] = eig(jac);
[Y,i] = min(abs(diag(D)));
V = V(:,i);
D = D(i,i);

if homds.nneg
    if homds.npos
        Q1U = Q0U * [eye(size(YU,1)); YU];
    else
        Q1U = [];
    end
    Qtot = [Q1U V];
    [Qn,Rn] = qr(Qtot);
    for i=1:homds.nneg
        f(end+1,1) = (ups(:,1) - x0)' * Qn(:,end-i+1);
    end
end
if homds.npos
    if homds.nneg
        Q1S = Q0S * [eye(size(YS,1)); YS];
    else
        Q1S = [];
    end
    Qtot = [Q1S V];
    [Qn,Rn] = qr(Qtot);
    for i=1:homds.npos
        f(end+1,1) = (ups(:,end) - x0)' * Qn(:,end-i+1);
    end
end

% Component 6
% -----------
% Distances from endpoints to equilibrium equal to epsilons
f(end+1,1) = norm(ups(:,1) - x0) - eps0;
f(end+1,1) = norm(ups(:,end) - x0) - eps1;
