function f = BVP_HTHom(x,x0,p,up,sp,T,eps1,YS,YU)

global HTHomds 

if HTHomds.index == 1
    % extract ups
    ups = reshape(x,HTHomds.nphase,HTHomds.tps);

    % Cycle-component
    % ---------------
    % Here the cycle is approximated in each mesh interval as a polynomial
    % of degree nCol, using the Lagrange polynomials. We approximate using the
    % collocation points in the subfunction below, and here we compute the
    % values of the polynomials at the mesh points.
    range1 = HTHomds.cols_p1;
    range2 = HTHomds.phases;
    f = [];
    p = HTHomds.P0;
    p = num2cell(p);
    for j=HTHomds.tsts
        % value of polynomial on each collocation point
        xval = ups(:,range1)*HTHomds.wt;  
        % derivative of polynomial on each collocation point
        derval  = ups(:,range1)*HTHomds.wpvec/HTHomds.dt(j);
        % evaluate function value on each collocation point
        for c=HTHomds.cols
            f(range2,1) = derval(:,c) - 2*HTHomds.T * feval(HTHomds.func, 0, xval(:,c), p{:});        
            % Shift for storage
            range2 = range2+HTHomds.nphase;
        end
        range1 = range1+HTHomds.ncol;
    end
    
    % Component 2
    % -----------
    % condition start_point
    Q0U = HTHomds.oldUnstableQ;
    for k = 1:HTHomds.nneg
        f(end+1,1) = Q0U(:,end-k+1)'*(ups(:,1)-HTHomds.x0);
    end

    % Component 3
    % -----------
    %condition SParams
    Q0S = HTHomds.oldStableQ;
    for k = 1:HTHomds.npos
        f(end+1,1) = sp(k) - 1/eps1*(ups(:,end)-HTHomds.x0)'*Q0S(:,HTHomds.nneg+k);
    end

    % Component 4
    % -----------
    % Distances from endpoint to equilibrium equal to epsilon1
    f(end+1,1) = norm(ups(:,1) - HTHomds.x0) - HTHomds.eps0;    
    f(end+1,1) = norm(ups(:,end) - HTHomds.x0) - eps1;

    % Component 5
    % -----------
    Q0U = HTHomds.oldUnstableQ;
    for k = 1:HTHomds.npos
        f(end+1,1) = HTHomds.eps0*up(k) - Q0U(:,k)'*(ups(:,1)-HTHomds.x0);
    end
    
elseif HTHomds.index == 2
    ups = reshape(x,HTHomds.nphase,HTHomds.tps);
    p = num2cell(p);

    % Cycle-component
    % ---------------
    % Here the cycle is approximated in each mesh interval as a polynomial
    % of degree nCol, using the Lagrange polynomials. We approximate using the
    % collocation points in the subfunction below, and here we compute the
    % values of the polynomials at the mesh points.
    range1 = HTHomds.cols_p1;
    range2 = HTHomds.phases;
    f = [];
    for j=HTHomds.tsts
        % value of polynomial on each collocation point
        xval = ups(:,range1)*HTHomds.wt;  
        % derivative of polynomial on each collocation point
        derval  = ups(:,range1)*HTHomds.wpvec/HTHomds.dt(j);
        % evaluate function value on each collocation point
        for c=HTHomds.cols
            f(range2,1) = derval(:,c) - 2*HTHomds.T * feval(HTHomds.func, 0, xval(:,c), p{:});        
            % Shift for storage
            range2 = range2+HTHomds.nphase;
        end
        range1 = range1+HTHomds.ncol;
    end

    % Component 2
    % -----------
    % equilibrium condition
    f(end+1:end+HTHomds.nphase,1) = feval(HTHomds.func, 0, x0, p{:}); 

    % Component 3
    % -----------
    % Ricatti equations
    if HTHomds.npos~=0 && HTHomds.nneg~=0
        f(end+1:end+2*HTHomds.Ysize) = ricattiEvalHTHom(x0,p,YU,YS);
    end        

    % Component 4
    % -----------
    % first vector along unstable eigenspaces   
    Q0U = HTHomds.oldUnstableQ;  
    if HTHomds.nneg
        Q1U = Q0U * [-YU'; eye(size(YU,1))];
        for i=1:HTHomds.nneg
            f(end+1,1) = (ups(:,1) - x0)' * Q1U(:,end-i+1);
        end
    end

    % Component 5
    % -----------
    Q0S = HTHomds.oldStableQ;
    if HTHomds.npos
        Q1S = Q0S * [-YS'; eye(size(YS,1))];
        for k = 1:HTHomds.npos  
            if k == HTHomds.ActiveSParams(1)
                f(end+1,1) = sp - 1/eps1*(ups(:,end)-x0)'*Q1S(:,k);
            else 
                f(end+1,1) = HTHomds.SParams(k) - 1/eps1*(ups(:,end)-x0)'*Q1S(:,k);
            end
        end
    end

    % Component 6
    % -----------
    % Distances from endpoints to equilibrium equal to epsilons
    f(end+1,1) = norm(ups(:,1) - x0) - HTHomds.eps0;
    f(end+1,1) = norm(ups(:,end) - x0) - eps1;

else
    
    % extract ups
    ups = reshape(x,HTHomds.nphase,HTHomds.tps);
    p = num2cell(p);

    % Cycle-component
    % ---------------
    % Here the cycle is approximated in each mesh interval as a polynomial
    % of degree nCol, using the Lagrange polynomials. We approximate using the
    % collocation points in the subfunction below, and here we compute the
    % values of the polynomials at the mesh points.
    range1 = HTHomds.cols_p1;
    range2 = HTHomds.phases;
    f = [];
    for j=HTHomds.tsts
        % value of polynomial on each collocation point
        xval = ups(:,range1)*HTHomds.wt;  
        % derivative of polynomial on each collocation point
        derval  = ups(:,range1)*HTHomds.wpvec/HTHomds.dt(j);
        % evaluate function value on each collocation point
        for c=HTHomds.cols
            f(range2,1) = derval(:,c) - 2*T * feval(HTHomds.func, 0, xval(:,c), p{:});        
            % Shift for storage
            range2 = range2+HTHomds.nphase;
        end
        range1 = range1+HTHomds.ncol;
    end

    % Component 2
    % -----------
    % equilibrium condition
    f(end+1:end+HTHomds.nphase,1) = feval(HTHomds.func, 0, x0, p{:}); 
    
    % Component 3
    % -----------
    % Ricatti equations
    if HTHomds.npos ~=0 && HTHomds.nneg ~= 0
        f(end+1:end+2*HTHomds.Ysize) = ricattiEvalHTHom(x0,p,YU,YS);
    end

    % Component 4
    % -----------
    % Last and first vectors along stable and unstable eigenspaces
    Q0S = HTHomds.oldStableQ;
    Q0U = HTHomds.oldUnstableQ;
    if HTHomds.nneg
        Q1U = Q0U * [-YU'; eye(size(YU,1))];
        for i=1:HTHomds.nneg
            f(end+1,1) = (ups(:,1) - x0)' * Q1U(:,end-i+1);
        end
    end
    if HTHomds.npos
        Q1S = Q0S * [-YS'; eye(size(YS,1))];
        for i=1:HTHomds.npos
            f(end+1,1) = (ups(:,end) - x0)' * Q1S(:,end-i+1);
        end
    end

    % Component 5
    % -----------
    % Distances from endpoints to equilibrium equal to epsilons
    f(end+1,1) = norm(ups(:,1) - x0) - HTHomds.eps0;
    f(end+1,1) = norm(ups(:,end) - x0) - eps1;

end

