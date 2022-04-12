function f = BVP_HTHet(x,x0,x1,p,up,sp,T,eps1,YS,YU)

%up bevat alle c's
%sp bevat alle tau's
global HTHetds 

if HTHetds.index == 1
    % extract ups
    ups = reshape(x,HTHetds.nphase,HTHetds.tps);

    % Cycle-component
    % ---------------
    % Here the cycle is approximated in each mesh interval as a polynomial
    % of degree nCol, using the Lagrange polynomials. We approximate using the
    % collocation points in the subfunction below, and here we compute the
    % values of the polynomials at the mesh points.
    range1 = HTHetds.cols_p1;
    range2 = HTHetds.phases;
    f = [];
    p = HTHetds.P0;
    p = num2cell(p);
    for j=HTHetds.tsts
        % value of polynomial on each collocation point
        xval = ups(:,range1)*HTHetds.wt;  
        % derivative of polynomial on each collocation point
        derval  = ups(:,range1)*HTHetds.wpvec/HTHetds.dt(j);
        % evaluate function value on each collocation point
        for c=HTHetds.cols
            f(range2,1) = derval(:,c) - 2*HTHetds.T * feval(HTHetds.func, 0, xval(:,c), p{:});        
            % Shift for storage
            range2 = range2+HTHetds.nphase;
        end
        range1 = range1+HTHetds.ncol;
    end
    
    % Component 2
    % -----------
    % condition start_point
    Q0U = HTHetds.oldUnstableQ;
    if (HTHetds.nphase-HTHetds.npos)
        for k = 1:HTHetds.nphase-HTHetds.npos
            f(end+1,1) = Q0U(:,end-k+1)'*(ups(:,1)-HTHetds.x0);
        end
    end

    % Component 3
    % -----------
    %condition SParams
    Q0S = HTHetds.oldStableQ;
    if (HTHetds.nphase-HTHetds.nneg)
        for k = 1:HTHetds.nphase-HTHetds.nneg
            f(end+1,1) = sp(k) - 1/eps1*(ups(:,end)-HTHetds.x1)'*Q0S(:,HTHetds.nneg+k);
        end
    end

    % Component 4
    % -----------
    % Distances from endpoint to equilibrium equal to epsilon1
    f(end+1,1) = norm(ups(:,1) - HTHetds.x0) - HTHetds.eps0;    
    f(end+1,1) = norm(ups(:,end) - HTHetds.x1) - eps1;

    % Component 5
    % -----------
    Q0U = HTHetds.oldUnstableQ;
    if HTHetds.npos
        for k = 1:HTHetds.npos
            f(end+1,1) = HTHetds.eps0*up(k) - Q0U(:,k)'*(ups(:,1)-HTHetds.x0);
        end
    end
    
elseif HTHetds.index == 2
    ups = reshape(x,HTHetds.nphase,HTHetds.tps);
    p = num2cell(p);

    % Cycle-component
    % ---------------
    % Here the cycle is approximated in each mesh interval as a polynomial
    % of degree nCol, using the Lagrange polynomials. We approximate using the
    % collocation points in the subfunction below, and here we compute the
    % values of the polynomials at the mesh points.
    range1 = HTHetds.cols_p1;
    range2 = HTHetds.phases;
    f = [];
    for j=HTHetds.tsts
        % value of polynomial on each collocation point
        xval = ups(:,range1)*HTHetds.wt;  
        % derivative of polynomial on each collocation point
        derval  = ups(:,range1)*HTHetds.wpvec/HTHetds.dt(j);
        % evaluate function value on each collocation point
        for c=HTHetds.cols
            f(range2,1) = derval(:,c) - 2*HTHetds.T * feval(HTHetds.func, 0, xval(:,c), p{:});        
            % Shift for storage
            range2 = range2+HTHetds.nphase;
        end
        range1 = range1+HTHetds.ncol;
    end

    % Component 2
    % -----------
    % equilibrium condition
    f(end+1:end+HTHetds.nphase,1) = feval(HTHetds.func, 0, x0, p{:}); 
    f(end+1:end+HTHetds.nphase,1) = feval(HTHetds.func, 0, x1, p{:});

    % Component 3
    % -----------
    % Ricatti equations
    if HTHetds.npos && (HTHetds.nphase-HTHetds.npos)
        f(end+1:end+(HTHetds.nphase-HTHetds.npos)*HTHetds.npos) = RicattiEvalHTHet(x0,p,1,YU);
    end
    if HTHetds.nneg && (HTHetds.nphase-HTHetds.nneg)
        f(end+1:end+(HTHetds.nphase-HTHetds.nneg)*HTHetds.nneg) = RicattiEvalHTHet(x1,p,0,YS);    
    end

    % Component 4
    % -----------
    % first vector along unstable eigenspaces   
    Q0U = HTHetds.oldUnstableQ;  
    if (HTHetds.nphase-HTHetds.npos)
        Q1U = Q0U * [-YU'; eye(size(YU,1))];
        for i=1:HTHetds.nphase-HTHetds.npos
            f(end+1,1) = (ups(:,1) - x0)' * Q1U(:,end-i+1);
        end
    end

    % Component 5
    % -----------
    Q0S = HTHetds.oldStableQ;
    if (HTHetds.nphase-HTHetds.nneg)
        Q1S = Q0S * [-YS'; eye(size(YS,1))];
        for k = 1:HTHetds.nphase-HTHetds.nneg                    
            f(end+1,1) = sp(k) - 1/eps1*(ups(:,end)-x1)'*Q1S(:,k);            
        end
    end

    % Component 6
    % -----------
    % Distances from endpoints to equilibrium equal to epsilons
    f(end+1,1) = norm(ups(:,1) - x0) - HTHetds.eps0;
    f(end+1,1) = norm(ups(:,end) - x1) - eps1;

else
    
    % extract ups
    ups = reshape(x,HTHetds.nphase,HTHetds.tps);
    p = num2cell(p);

    % Cycle-component
    % ---------------
    % Here the cycle is approximated in each mesh interval as a polynomial
    % of degree nCol, using the Lagrange polynomials. We approximate using the
    % collocation points in the subfunction below, and here we compute the
    % values of the polynomials at the mesh points.
    range1 = HTHetds.cols_p1;
    range2 = HTHetds.phases;
    f = [];
    for j=HTHetds.tsts
        % value of polynomial on each collocation point
        xval = ups(:,range1)*HTHetds.wt;  
        % derivative of polynomial on each collocation point
        derval  = ups(:,range1)*HTHetds.wpvec/HTHetds.dt(j);
        % evaluate function value on each collocation point
        for c=HTHetds.cols
            f(range2,1) = derval(:,c) - 2*T * feval(HTHetds.func, 0, xval(:,c), p{:});        
            % Shift for storage
            range2 = range2+HTHetds.nphase;
        end
        range1 = range1+HTHetds.ncol;
    end

    % Component 2
    % -----------
    % equilibrium condition
    f(end+1:end+HTHetds.nphase,1) = feval(HTHetds.func, 0, x0, p{:}); 
    f(end+1:end+HTHetds.nphase,1) = feval(HTHetds.func, 0, x1, p{:}); 
    
    % Component 3
    % -----------
    % Ricatti equations
    if HTHetds.npos && (HTHetds.nphase-HTHetds.npos)
        f(end+1:end+(HTHetds.nphase-HTHetds.npos)*HTHetds.npos) = RicattiEvalHTHet(x0,p,1,YU);
    end
    if HTHetds.nneg && (HTHetds.nphase-HTHetds.nneg)
        f(end+1:end+(HTHetds.nphase-HTHetds.nneg)*HTHetds.nneg) = RicattiEvalHTHet(x1,p,0,YS);    
    end
    
    % Component 4
    % -----------
    % Last and first vectors along stable and unstable eigenspaces
    Q0S = HTHetds.oldStableQ;
    Q0U = HTHetds.oldUnstableQ;
    if (HTHetds.nphase-HTHetds.npos)
        Q1U = Q0U * [-YU'; eye(size(YU,1))];
        for i=1:HTHetds.nphase-HTHetds.npos
            f(end+1,1) = (ups(:,1) - x0)' * Q1U(:,end-i+1);
        end
    end
    if (HTHetds.nphase-HTHetds.nneg)
        Q1S = Q0S * [-YS'; eye(size(YS,1))];
        for i=1:HTHetds.nphase-HTHetds.nneg
            f(end+1,1) = (ups(:,end) - x1)' * Q1S(:,end-i+1);
        end
    end

    % Component 5
    % -----------
    % Distances from endpoints to equilibrium equal to epsilons
    f(end+1,1) = norm(ups(:,1) - x0) - HTHetds.eps0;
    f(end+1,1) = norm(ups(:,end) - x1) - eps1;

end

