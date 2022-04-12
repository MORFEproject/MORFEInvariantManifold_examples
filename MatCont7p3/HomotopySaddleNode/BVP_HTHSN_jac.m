% Build_connec2ds up jacobian of boundary value problem
%
% ============================================

function result = BVP_HTHSN_jac(odefile1,x,up,sp,T,eps1)
global HTHSNds

if HTHSNds.index == 1

    ups = reshape(x,HTHSNds.nphase,HTHSNds.tps);
    tmpperiod = 2 * HTHSNds.T;

    % Allocate space for sparse jacobian
    result = spalloc(HTHSNds.ncoords+HTHSNds.npos+2,HTHSNds.ncoords+HTHSNds.npos+3,...
    (HTHSNds.ncoords+1)^2);
    range1 = HTHSNds.cols_p1;
    range2 = 1:HTHSNds.ncol*HTHSNds.nphase;
    range3 = 1:(HTHSNds.ncol+1)*HTHSNds.nphase;

    % Component 1 
    % ===========
    p = HTHSNds.P0;
    p = num2cell(p);
    for i=HTHSNds.tsts
        range4 = HTHSNds.phases;
        % value of polynomial in each mesh point
        xval = ups(:,range1) * HTHSNds.wt;
        wploc = HTHSNds.wp/HTHSNds.dt(i);
        % evaluate part of Jacobian matrix
        for j=HTHSNds.cols
            xtmp = xval(:,j);  
        
            jac=cjac(HTHSNds.func,HTHSNds.Jacobian,xtmp,p,[]);
            sysjac(range4,:) = fastkron(HTHSNds.ncol,HTHSNds.nphase,HTHSNds.wt(:,j)',jac);
        
            range4 = range4+HTHSNds.nphase;
        end
        % Store result

        result(range2,range3) = [wploc-tmpperiod*sysjac];

        % shift to following intervals
        range1 = range1 + HTHSNds.ncol;
        range2 = range2 + HTHSNds.ncol*HTHSNds.nphase;
        range3 = range3 + HTHSNds.ncol*HTHSNds.nphase;
    end

    % Component 2 
    % ===========
    Q0U = HTHSNds.oldUnstableQ;
    jac = cjac(HTHSNds.func,HTHSNds.Jacobian,HTHSNds.x0,p,[]);
    [V,D] = eig(jac);
    [Y,i] = min(abs(diag(D)));
    V = V(:,i);
    Qtot = [Q0U(:,1:HTHSNds.npos) V];
    [Qn,Rn] = qr(Qtot);
    if HTHSNds.nneg
        for i=1:HTHSNds.nneg
            result(HTHSNds.ncoords-HTHSNds.nphase+i,1:HTHSNds.nphase) = Qn(:,end-i+1)';    
        end
    end

    % Component 3
    % ===========
    Q0S = HTHSNds.oldStableQ;
    Qtot = [Q0S(:,1:HTHSNds.nneg) V];
    [Qn,Rn] = qr(Qtot);
    if HTHSNds.npos
        for i = 1:HTHSNds.npos
            result(HTHSNds.ncoords-HTHSNds.nphase+HTHSNds.nneg+i,HTHSNds.ncoords-HTHSNds.nphase+1:HTHSNds.ncoords) = -Qn(:,HTHSNds.nneg+1+i)';       
            result(HTHSNds.ncoords-HTHSNds.nphase+HTHSNds.nneg+i,end) = 1/eps1^2*(ups(:,end)-HTHSNds.x0)'*Qn(:,HTHSNds.nneg+1+i);       
        end
    end

    for i = 1:size(HTHSNds.ActiveSParams,2)
        index = HTHSNds.ActiveSParams(i);
        result(HTHSNds.ncoords-HTHSNds.nphase+HTHSNds.nneg+index,HTHSNds.ncoords+size(HTHSNds.ActiveUParams,2)+i) = 1;
    end

    % Component 4
    % ===========
    val = ups(:,1) - HTHSNds.x0;
    val = val' / norm(val);
    result(HTHSNds.ncoords,1:HTHSNds.nphase) = val;

    val = ups(:,end) - HTHSNds.x0;
    val = val' / norm(val);
    result(HTHSNds.ncoords+1,HTHSNds.ncoords-HTHSNds.nphase+1:HTHSNds.ncoords) = val;   
    result(HTHSNds.ncoords+1,end) = -1;    

    % Component 5
    % ===========
    for k = 1:size(HTHSNds.ActiveUParams,2)
        index = HTHSNds.ActiveUParams(k);
        result(HTHSNds.ncoords+1+index,HTHSNds.ncoords+k) = HTHSNds.eps0;      
    end

    Qtot = [Q0U(:,1:HTHSNds.npos) V];
    [Qn,Rn] = qr(Qtot);
    for k = 1:HTHSNds.npos+1
        result(HTHSNds.ncoords+1+k,1:HTHSNds.nphase) = -Qn(:,k)';      
    end

else 
    ups = reshape(x,HTHSNds.nphase,HTHSNds.tps);
    p = HTHSNds.P0;
    p = num2cell(p);

    pars = HTHSNds.ncoords +1;
    tmpperiod = 2 * T;

    % Allocate space for sparse jacobian
    result = spalloc(HTHSNds.ncoords+1,...
        HTHSNds.ncoords+2,(HTHSNds.ncoords+1)^2);
    range1 = HTHSNds.cols_p1;
    range2 = 1:HTHSNds.ncol*HTHSNds.nphase;
    range3 = 1:(HTHSNds.ncol+1)*HTHSNds.nphase;

    % Component 1 
    % ===========
    for i=HTHSNds.tsts
        range4 = HTHSNds.phases;
        % value of polynomial in each mesh point
        xval = ups(:,range1) * HTHSNds.wt;
        wploc = HTHSNds.wp/HTHSNds.dt(i);
        % evaluate part of Jacobian matrix
        for j=HTHSNds.cols
            xtmp = xval(:,j);
            frhs(range4) = feval(odefile1, 0, xtmp, p{:});
            jac=cjac(HTHSNds.func,HTHSNds.Jacobian,xtmp,p,[]);
            sysjac(range4,:) = fastkron(HTHSNds.ncol,HTHSNds.nphase,HTHSNds.wt(:,j)',jac);
        
            range4 = range4+HTHSNds.nphase;
        end
        % Store result
        result(range2,[range3 pars]) = [wploc-tmpperiod*sysjac   -2*frhs'];

        % shift to following intervals
        range1 = range1 + HTHSNds.ncol;
        range2 = range2 + HTHSNds.ncol*HTHSNds.nphase;
        range3 = range3 + HTHSNds.ncol*HTHSNds.nphase;
    end
      
    % Component 2 (Last and first vectors along stable and unstable eigenspaces)
    % ===========
    Q0U = HTHSNds.oldUnstableQ;  
    jac = cjac(HTHSNds.func,HTHSNds.Jacobian,HTHSNds.x0,p,[]);
    [V,D] = eig(jac);
    [Y,i] = min(abs(diag(D)));
    V = V(:,i);
   
    if HTHSNds.nneg
        Qtot = [Q0U(:,1:HTHSNds.npos) V];
        [Qn,Rn] = qr(Qtot);
        for i = 1:HTHSNds.nneg
            result(HTHSNds.ncoords-HTHSNds.nphase+i,1:HTHSNds.nphase) = Qn(:,end-i+1)';       
        end
    end

    Q0S = HTHSNds.oldStableQ;
       
    if HTHSNds.npos
        Qtot = [Q0S(:,1:HTHSNds.nneg) V];
        [Qn,Rn] = qr(Qtot);
        for i = 1:HTHSNds.npos
            result(HTHSNds.ncoords-HTHSNds.nphase+HTHSNds.nneg+i,HTHSNds.ncoords-HTHSNds.nphase+1:HTHSNds.ncoords) = Qn(:,end-i+1)';       
        end
    end

    % Component 6 (Distances from endpoints to equilibrium equal to epsilons)
    % ===========
    val = ups(:,1) - HTHSNds.x0;
    val = val' / norm(val);
    result(end-1,HTHSNds.phases) = val;
   
    val = ups(:,end) - HTHSNds.x0;
    val = val' / norm(val);
    result(end,HTHSNds.ncoords-HTHSNds.nphase+1:HTHSNds.ncoords) = val;
    result(end,HTHSNds.ncoords+2) = -1;

end