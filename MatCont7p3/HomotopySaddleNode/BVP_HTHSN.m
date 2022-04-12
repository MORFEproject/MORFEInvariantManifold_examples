function f = BVP_HTHSN(x,up,sp,T,eps1)

global HTHSNds 

if HTHSNds.index == 1
    % extract ups
    ups = reshape(x,HTHSNds.nphase,HTHSNds.tps);

    % Cycle-component
    % ---------------
    % Here the cycle is approximated in each mesh interval as a polynomial
    % of degree nCol, using the Lagrange polynomials. We approximate using the
    % collocation points in the subfunction below, and here we compute the
    % values of the polynomials at the mesh points.
    range1 = HTHSNds.cols_p1;
    range2 = HTHSNds.phases;
    f = [];
    p = HTHSNds.P0;
    p = num2cell(p);
    for j=HTHSNds.tsts
        % value of polynomial on each collocation point
        xval = ups(:,range1)*HTHSNds.wt;  
        % derivative of polynomial on each collocation point
        derval  = ups(:,range1)*HTHSNds.wpvec/HTHSNds.dt(j);
        % evaluate function value on each collocation point
        for c=HTHSNds.cols
            f(range2,1) = derval(:,c) - 2*HTHSNds.T * feval(HTHSNds.func, 0, xval(:,c), p{:});        
            % Shift for storage
            range2 = range2+HTHSNds.nphase;
        end
        range1 = range1+HTHSNds.ncol;
    end
    
    % Component 2
    % -----------
    % condition start_point
    Q0U = HTHSNds.oldUnstableQ;
    jac = cjac(HTHSNds.func,HTHSNds.Jacobian,HTHSNds.x0,p,[]);
    [V,D] = eig(jac);
    [Y,i] = min(abs(diag(D)));
    V = V(:,i);
    Qtot = [Q0U(:,1:HTHSNds.npos) V];
    [Qn,Rn] = qr(Qtot);
    
    for k = 1:HTHSNds.nneg
        f(end+1,1) = Qn(:,end-k+1)'*(ups(:,1)-HTHSNds.x0);
    end

    % Component 3
    % -----------
    %condition SParams
    Q0S = HTHSNds.oldStableQ;
    Qtot = [Q0S(:,1:HTHSNds.nneg) V];
    [Qn,Rn] = qr(Qtot);
    
    for k = 1:HTHSNds.npos
        f(end+1,1) = sp(k) - 1/eps1*(ups(:,end)-HTHSNds.x0)'*Qn(:,HTHSNds.nneg+1+k);
    end

    % Component 4
    % -----------
    % Distances from endpoint to equilibrium equal to epsilon1
    f(end+1,1) = norm(ups(:,1) - HTHSNds.x0) - HTHSNds.eps0;    
    f(end+1,1) = norm(ups(:,end) - HTHSNds.x0) - eps1;

    % Component 5
    % -----------
    Q0U = HTHSNds.oldUnstableQ;
    Qtot = [Q0U(:,1:HTHSNds.npos) V];
    [Qn,Rn] = qr(Qtot);
    
    for k = 1:HTHSNds.npos+1
        f(end+1,1) = HTHSNds.eps0*up(k) - Qn(:,k)'*(ups(:,1)-HTHSNds.x0);
    end

else
    
    % extract ups
    ups = reshape(x,HTHSNds.nphase,HTHSNds.tps);
    p = HTHSNds.P0;
    p = num2cell(p);

    % Cycle-component
    % ---------------
    % Here the cycle is approximated in each mesh interval as a polynomial
    % of degree nCol, using the Lagrange polynomials. We approximate using the
    % collocation points in the subfunction below, and here we compute the
    % values of the polynomials at the mesh points.
    range1 = HTHSNds.cols_p1;
    range2 = HTHSNds.phases;
    f = [];
    for j=HTHSNds.tsts
        % value of polynomial on each collocation point
        xval = ups(:,range1)*HTHSNds.wt;  
        % derivative of polynomial on each collocation point
        derval  = ups(:,range1)*HTHSNds.wpvec/HTHSNds.dt(j);
        % evaluate function value on each collocation point
        for c=HTHSNds.cols
            f(range2,1) = derval(:,c) - 2*T * feval(HTHSNds.func, 0, xval(:,c), p{:});        
            % Shift for storage
            range2 = range2+HTHSNds.nphase;
        end
        range1 = range1+HTHSNds.ncol;
    end

    % Component 2
    % -----------
    % Last and first vectors along stable and unstable eigenspaces
    Q0S = HTHSNds.oldStableQ;
    Q0U = HTHSNds.oldUnstableQ;
    jac = cjac(HTHSNds.func,HTHSNds.Jacobian,HTHSNds.x0,p,[]);
    [V,D] = eig(jac);
    [Y,i] = min(abs(diag(D)));
    V = V(:,i);
       
    if HTHSNds.nneg
        Qtot = [Q0U(:,1:HTHSNds.npos) V];
        [Qn,Rn] = qr(Qtot);
        for i=1:HTHSNds.nneg
            f(end+1,1) = (ups(:,1) - HTHSNds.x0)' * Qn(:,end-i+1);
        end
    end
    
    if HTHSNds.npos        
        Qtot = [Q0S(:,1:HTHSNds.nneg) V];
        [Qn,Rn] = qr(Qtot);
        for i=1:HTHSNds.npos
            f(end+1,1) = (ups(:,end) - HTHSNds.x0)' * Qn(:,end-i+1);
        end
    end

    % Component 5
    % -----------
    % Distances from endpoints to equilibrium equal to epsilons
    f(end+1,1) = norm(ups(:,1) - HTHSNds.x0) - HTHSNds.eps0;
    f(end+1,1) = norm(ups(:,end) - HTHSNds.x0) - eps1;

end

