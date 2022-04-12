% Build_connec2ds up jacobian of boundary value problem
%
% ============================================

function result = BVP_HTHet_jac(odefile1,x,x0,x1,p,up,sp,T,eps1,YS,YU)
global HTHetds

if HTHetds.index == 1
    ups = reshape(x,HTHetds.nphase,HTHetds.tps);
    tmpperiod = 2 * HTHetds.T;

    % Allocate space for sparse jacobian
    result = spalloc(HTHetds.ncoords+HTHetds.nphase-HTHetds.nneg+2,HTHetds.ncoords+HTHetds.nphase-HTHetds.nneg+3,...
    (HTHetds.ncoords+1)^2);
    range1 = HTHetds.cols_p1;
    range2 = 1:HTHetds.ncol*HTHetds.nphase;
    range3 = 1:(HTHetds.ncol+1)*HTHetds.nphase;

    % Component 1 
    % ===========
    p = HTHetds.P0;
    p = num2cell(p);
    for i=HTHetds.tsts
        range4 = HTHetds.phases;
        % value of polynomial in each mesh point
        xval = ups(:,range1) * HTHetds.wt;
        wploc = HTHetds.wp/HTHetds.dt(i);
        % evaluate part of Jacobian matrix
        for j=HTHetds.cols
            xtmp = xval(:,j);  
        
            jac=cjac(HTHetds.func,HTHetds.Jacobian,xtmp,p,[]);
            sysjac(range4,:) = fastkron(HTHetds.ncol,HTHetds.nphase,HTHetds.wt(:,j)',jac);
        
            range4 = range4+HTHetds.nphase;
        end
        % Store result

        result(range2,range3) = [wploc-tmpperiod*sysjac];

        % shift to following intervals
        range1 = range1 + HTHetds.ncol;
        range2 = range2 + HTHetds.ncol*HTHetds.nphase;
        range3 = range3 + HTHetds.ncol*HTHetds.nphase;
    end

    % Component 2 
    % ===========
    Q0U = HTHetds.oldUnstableQ;
    if (HTHetds.nphase-HTHetds.npos)
        for i=1:HTHetds.nphase-HTHetds.npos
            result(HTHetds.ncoords-HTHetds.nphase+i,1:HTHetds.nphase) = Q0U(:,end-i+1)';    
        end
    end

    % Component 3
    % ===========
    QS = HTHetds.oldStableQ;
    if HTHetds.nphase-HTHetds.nneg
        for i = 1:HTHetds.nphase-HTHetds.nneg
            result(HTHetds.ncoords-HTHetds.npos+i,HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords) = -1/eps1*QS(:,HTHetds.nneg+i)';       
            result(HTHetds.ncoords-HTHetds.npos+i,end) = 1/eps1^2*(ups(:,end)-HTHetds.x1)'*QS(:,HTHetds.nneg+i);       
        end
    end

    for i = 1:size(HTHetds.ActiveSParams,2)
        index = HTHetds.ActiveSParams(i);
        result(HTHetds.ncoords-HTHetds.npos+index,HTHetds.ncoords+size(HTHetds.ActiveUParams,2)+i) = 1;
    end

    % Component 4
    % ===========
    val = ups(:,1) - HTHetds.x0;
    val = val' / norm(val);
    result(HTHetds.ncoords+HTHetds.nphase-HTHetds.npos-HTHetds.nneg+1,1:HTHetds.nphase) = val;

    val = ups(:,end) - HTHetds.x1;
    val = val' / norm(val);
    result(HTHetds.ncoords+HTHetds.nphase-HTHetds.npos-HTHetds.nneg+2,HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords) = val;
    
    result(HTHetds.ncoords+HTHetds.nphase-HTHetds.npos-HTHetds.nneg+2,end) = -1;    

    % Component 5
    % ===========
    for k = 1:size(HTHetds.ActiveUParams,2)
        index = HTHetds.ActiveUParams(k);
        result(HTHetds.ncoords+HTHetds.nphase-HTHetds.npos-HTHetds.nneg+2+index,HTHetds.ncoords+k) = HTHetds.eps0;      
    end

    if HTHetds.npos
        for k = 1:HTHetds.npos      
            result(HTHetds.ncoords+HTHetds.nphase-HTHetds.npos-HTHetds.nneg+2+k,1:HTHetds.nphase) = -Q0U(:,k)';      
        end
    end

elseif HTHetds.index == 2
   
    ups = reshape(x,HTHetds.nphase,HTHetds.tps);
    p = num2cell(p);
    
    pars = HTHetds.PeriodIdx + (1:length(HTHetds.ActiveParams));
    tmpperiod = 2 * HTHetds.T;

    % Allocate space for sparse jacobian
    result = spalloc(HTHetds.ncoords+HTHetds.YUsize+HTHetds.YSsize+3*HTHetds.nphase-HTHetds.npos-HTHetds.nneg+2,...
        HTHetds.ncoords+3*HTHetds.nphase-HTHetds.npos-HTHetds.nneg+3+HTHetds.YUsize+HTHetds.YSsize,(HTHetds.ncoords+1)^2);
    range1 = HTHetds.cols_p1;
    range2 = 1:HTHetds.ncol*HTHetds.nphase;
    range3 = 1:(HTHetds.ncol+1)*HTHetds.nphase;

    % Component 1 
    % ===========
    for i=HTHetds.tsts
        range4 = HTHetds.phases;
        % value of polynomial in each mesh point
        xval = ups(:,range1) * HTHetds.wt;
        wploc = HTHetds.wp/HTHetds.dt(i);
        % evaluate part of Jacobian matrix
        for j=HTHetds.cols
            xtmp = xval(:,j);
            jac=cjac(HTHetds.func,HTHetds.Jacobian,xtmp,p,HTHetds.ActiveParams);
            sysjac(range4,:) = fastkron(HTHetds.ncol,HTHetds.nphase,HTHetds.wt(:,j)',jac);
            sysjacp(range4,:) = cjacp(HTHetds.func,HTHetds.JacobianP,xtmp,p,HTHetds.ActiveParams);
        
            range4 = range4+HTHetds.nphase;
        end
        % Store result    
        result(range2,[range3 pars]) = [wploc-tmpperiod*sysjac    -tmpperiod*sysjacp];
    
        % shift to following intervals
        range1 = range1 + HTHetds.ncol;
        range2 = range2 + HTHetds.ncol*HTHetds.nphase;
        range3 = range3 + HTHetds.ncol*HTHetds.nphase;
    end

    % Component 2 (equilibrium)
    % ===========
    result(HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords, HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase) = cjac(HTHetds.func,HTHetds.Jacobian,x0,p,HTHetds.ActiveParams);
    result(HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords, HTHetds.PeriodIdx+1:HTHetds.PeriodIdx+length(HTHetds.ActiveParams)) = cjacp(HTHetds.func,HTHetds.JacobianP,x0,p,HTHetds.ActiveParams);
    result(HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase, HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+2*HTHetds.nphase) = cjac(HTHetds.func,HTHetds.Jacobian,x1,p,HTHetds.ActiveParams);
    result(HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase, HTHetds.PeriodIdx+1:HTHetds.PeriodIdx+length(HTHetds.ActiveParams)) = cjacp(HTHetds.func,HTHetds.JacobianP,x1,p,HTHetds.ActiveParams);    

    % Component 3 (Eigenspaces)
    % ===========
    Q0U = HTHetds.oldUnstableQ;
    Q0S = HTHetds.oldStableQ;
    hess0 = chess(HTHetds.func,HTHetds.Jacobian,HTHetds.Hessians,x0,p,HTHetds.ActiveParams);
    hessp0 = chessp(HTHetds.func,HTHetds.Jacobian,HTHetds.HessiansP,x0,p,HTHetds.ActiveParams);
    hess1 = chess(HTHetds.func,HTHetds.Jacobian,HTHetds.Hessians,x1,p,HTHetds.ActiveParams);
    hessp1 = chessp(HTHetds.func,HTHetds.Jacobian,HTHetds.HessiansP,x1,p,HTHetds.ActiveParams);
    if HTHetds.YUsize ~= 0
    for j=1:HTHetds.nphase+length(HTHetds.ActiveParams)
        if j<=HTHetds.nphase
            Atemp = hess0(:,:,j);
        else
            Atemp = hessp0(:,:,j-HTHetds.nphase);
        end 
        [U11, U12, UE21, U22] = RicattiCoeffHTHet(Q0U,Atemp,HTHetds.npos);
        tmp = (U22*YU - YU*U11 + UE21 - YU*U12*YU)';
        if j<=HTHetds.nphase
            result(HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize,HTHetds.ncoords+j) = tmp(:);
        else
            k = j-HTHetds.nphase;
            result(HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize,HTHetds.PeriodIdx+k) = tmp(:);
        end
    end
    end
    if HTHetds.YSsize ~= 0
    for j=1:HTHetds.nphase+length(HTHetds.ActiveParams)
        if j<=HTHetds.nphase
            Atemp = hess1(:,:,j);
        else
            Atemp = hessp1(:,:,j-HTHetds.nphase);
        end 
        [S11, S12, SE21, S22] = RicattiCoeffHTHet(Q0S,Atemp,HTHetds.nneg);
        tmp = (S22*YS - YS*S11 + SE21 - YS*S12*YS)';

        result(HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+1:HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+HTHetds.YSsize,HTHetds.ncoords+HTHetds.nphase+j) = tmp(:);        
    end
    end                           

    % UNSTABLE
    % ricatti blocks from unstable eigenspace
    if HTHetds.YUsize    
    A=cjac(HTHetds.func,HTHetds.Jacobian,x0,p,HTHetds.ActiveParams);
    derU = [];
    [U11, U12, UE21, U22] = RicattiCoeffHTHet(Q0U,A,HTHetds.npos);
    U12Y = U12 * YU;
    YU12 = YU * U12;
    for j=1:HTHetds.npos
        for i=1:HTHetds.nphase-HTHetds.npos
            myres = sparse(HTHetds.nphase-HTHetds.npos,HTHetds.npos,0);
            myres(:,j) = myres(:,j) + U22(:,i);
            myres(i,:) = myres(i,:) - U11(j,:);        
            myres(i,:) = myres(i,:) - U12Y(j,:);
            myres(:,j) = myres(:,j) - YU12(:,i);        
            % Derivative is computed to YU(i,j) i=1:HTHetds.nneg j=1:HTHetds.npos
            % must come at column (j-1)*(HTHetds.nneg)+i, 
            % in form T and col-first == row-first
            rows = size(myres,1);
            cols = size(myres,2);
            for k=1:rows
                derU((k-1)*cols+1:k*cols,(j-1)*(HTHetds.nphase-HTHetds.npos)+i) = myres(k,:)';
            end
        end
    end
    end

    % STABLE
    if HTHetds.YSsize
    A=cjac(HTHetds.func,HTHetds.Jacobian,x1,p,HTHetds.ActiveParams);
    derS = [];
    [S11, S12, SE21, S22] = RicattiCoeffHTHet(Q0S,A,HTHetds.nneg);
    S12Y = S12 * YS;
    YS12 = YS * S12;
    for i=1:HTHetds.nphase-HTHetds.nneg
        for j=1:HTHetds.nneg
            myres = sparse(HTHetds.nphase-HTHetds.nneg,HTHetds.nneg,0);
            myres(:,j) = myres(:,j) + S22(:,i);
            myres(i,:) = myres(i,:) - S11(j,:);
            myres(i,:) = myres(i,:) - S12Y(j,:);
            myres(:,j) = myres(:,j) - YS12(:,i);        
            % Derivative is computed to YS(i,j)
            % must come at column (j-1)*(HTHetds.npos)+i, 
            % in form T and col-first == row-first
            rows = size(myres,1);
            cols = size(myres,2);
            for k=1:rows
                derS((k-1)*cols+1:k*cols,(j-1)*(HTHetds.nphase-HTHetds.nneg)+i) = myres(k,:)';
            end
        end
    end
    end
    if HTHetds.YUsize
    result(HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize,end-HTHetds.YUsize-HTHetds.YSsize+1:end-HTHetds.YSsize) = derU;
    end
    if HTHetds.YSsize    
    result(HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+1:HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+HTHetds.YSsize,end-HTHetds.YSsize+1:end) = derS;
    end

    % Component 4 (Last and first vectors along stable and unstable eigenspaces)
    % ===========
    if (HTHetds.nphase-HTHetds.npos)
        indxs = [1:HTHetds.nphase   HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase];
        vect = ups(:,1) - x0;
        Q1U = Q0U * [-YU'; eye(size(YU,1))];
        vQ = -vect' * Q0U;
        for i=1:HTHetds.nphase-HTHetds.npos
            result(HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+HTHetds.YSsize+i,indxs) = [Q1U(:,end-i+1)'     -Q1U(:,end-i+1)'];
        end
        if HTHetds.npos
            myeye = zeros(HTHetds.nphase-HTHetds.npos);
            for j=1:HTHetds.nphase-HTHetds.npos
                myeye(j,HTHetds.nphase-HTHetds.npos-j+1) = 1;
            end
            for j=1:HTHetds.npos
                result(HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+HTHetds.YSsize+1:HTHetds.ncoords+2*HTHetds.nphase+HTHetds.YUsize+HTHetds.YSsize-HTHetds.npos,end-HTHetds.YUsize-HTHetds.YSsize+(j-1)*(HTHetds.nphase-HTHetds.npos)+1:end-HTHetds.YUsize-HTHetds.YSsize+j*(HTHetds.nphase-HTHetds.npos)) = vQ(:,j) * myeye;
            end
        end
    end

    % Component 5 
    % ===========
    Q0S = HTHetds.oldStableQ;
    if (HTHetds.nphase-HTHetds.nneg)
        Q1S = Q0S * [-YS'; eye(size(YS,1))];
        for i = 1:HTHetds.nphase-HTHetds.nneg
            result(end-2-HTHetds.nphase+HTHetds.nneg+i,HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords) = -1/eps1*Q1S(:,i)';    
            result(end-2-HTHetds.nphase+HTHetds.nneg+i,HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+2*HTHetds.nphase) = 1/eps1*Q1S(:,i)';       
        end
        if HTHetds.nneg
            vect = ups(:,end) - x1;
            vQ = vect' * Q0S;  
            myeye = zeros(HTHetds.nphase-HTHetds.nneg);
            for j=1:HTHetds.nphase-HTHetds.nneg
                myeye(j,j) = 1;
            end
            %afgeleide naar de elementen van YS
            for j=1:HTHetds.nneg
                result(end-2-HTHetds.nphase+HTHetds.nneg+1:end-2,end-HTHetds.YSsize+(j-1)*(HTHetds.nphase-HTHetds.nneg)+1:end-HTHetds.YSsize+j*(HTHetds.nphase-HTHetds.nneg)) = vQ(:,j) * myeye;
            end
        end

        vec = 1/eps1^2*(ups(:,end)-x1)'*Q1S;
        result(end-2-HTHetds.nphase+HTHetds.nneg+1:end-2,end-HTHetds.YSsize-HTHetds.YUsize) = vec';
    end
    for i = 1:length(HTHetds.ActiveSParams)
        k = HTHetds.ActiveSParams(i);
        result(end-2-HTHetds.nphase+HTHetds.nneg+k,HTHetds.PeriodIdx+length(HTHetds.ActiveParams)+i) = 1; %afgeleide naar actieve tau
    end

    % Component 6 (Distances from endpoints to equilibrium equal to epsilons)
    % ===========
    indxs = [HTHetds.phases   HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase];
    val = ups(:,1) - x0;
    val = val' / norm(val);
    result(end-1,indxs) = [val   -val];
    
    indxs = [HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords HTHetds.ncoords+HTHetds.nphase+1:HTHetds.PeriodIdx];
    val = ups(:,end) - x1;
    val = val' / norm(val);
    result(end,indxs) = [val   -val];

    result(end,end-HTHetds.YUsize-HTHetds.YSsize) = -1;

else 
    
    ups = reshape(x,HTHetds.nphase,HTHetds.tps);
    p = num2cell(p);

    pars = (0:length(HTHetds.ActiveParams))+HTHetds.PeriodIdx +1;
    tmpperiod = 2 * T;

    % Allocate space for sparse jacobian
    result = spalloc(HTHetds.ncoords+3*HTHetds.nphase-HTHetds.nneg-HTHetds.npos+2+HTHetds.YUsize+HTHetds.YSsize,...
        HTHetds.ncoords+3*HTHetds.nphase-HTHetds.nneg-HTHetds.npos+3+HTHetds.YUsize+HTHetds.YSsize,(HTHetds.ncoords+1)^2);
    range1 = HTHetds.cols_p1;
    range2 = 1:HTHetds.ncol*HTHetds.nphase;
    range3 = 1:(HTHetds.ncol+1)*HTHetds.nphase;

    % Component 1 
    % ===========
    for i=HTHetds.tsts
        range4 = HTHetds.phases;
        % value of polynomial in each mesh point
        xval = ups(:,range1) * HTHetds.wt;
        wploc = HTHetds.wp/HTHetds.dt(i);
        % evaluate part of Jacobian matrix
        for j=HTHetds.cols
            xtmp = xval(:,j);
            frhs(range4) = feval(odefile1, 0, xtmp, p{:});
            jac=cjac(HTHetds.func,HTHetds.Jacobian,xtmp,p,HTHetds.ActiveParams);
            sysjac(range4,:) = fastkron(HTHetds.ncol,HTHetds.nphase,HTHetds.wt(:,j)',jac);
            sysjacp(range4,:) = cjacp(HTHetds.func,HTHetds.JacobianP,xtmp,p,HTHetds.ActiveParams);
        
            range4 = range4+HTHetds.nphase;
        end
        % Store result
        result(range2,[range3 pars]) = [wploc-tmpperiod*sysjac    -tmpperiod*sysjacp    -2*frhs'];

        % shift to following intervals
        range1 = range1 + HTHetds.ncol;
        range2 = range2 + HTHetds.ncol*HTHetds.nphase;
        range3 = range3 + HTHetds.ncol*HTHetds.nphase;
    end

    % Component 2 (equilibrium)
    % ===========
    result(HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords, HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase) = cjac(HTHetds.func,HTHetds.Jacobian,x0,p,HTHetds.ActiveParams);
    result(HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords, HTHetds.PeriodIdx+1:HTHetds.PeriodIdx+length(HTHetds.ActiveParams)) = cjacp(HTHetds.func,HTHetds.JacobianP,x0,p,HTHetds.ActiveParams);
    result(HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase, HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+2*HTHetds.nphase) = cjac(HTHetds.func,HTHetds.Jacobian,x1,p,HTHetds.ActiveParams);
    result(HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase, HTHetds.PeriodIdx+1:HTHetds.PeriodIdx+length(HTHetds.ActiveParams)) = cjacp(HTHetds.func,HTHetds.JacobianP,x1,p,HTHetds.ActiveParams);    

    % Component 3 (Eigenspaces)
    % ===========
    Q0U = HTHetds.oldUnstableQ;
    Q0S = HTHetds.oldStableQ;
    hess0 = chess(HTHetds.func,HTHetds.Jacobian,HTHetds.Hessians,x0,p,HTHetds.ActiveParams);
    hessp0 = chessp(HTHetds.func,HTHetds.Jacobian,HTHetds.HessiansP,x0,p,HTHetds.ActiveParams);
    hess1 = chess(HTHetds.func,HTHetds.Jacobian,HTHetds.Hessians,x1,p,HTHetds.ActiveParams);
    hessp1 = chessp(HTHetds.func,HTHetds.Jacobian,HTHetds.HessiansP,x1,p,HTHetds.ActiveParams);
    if HTHetds.YUsize
    for j=1:HTHetds.nphase+length(HTHetds.ActiveParams)
        if j<=HTHetds.nphase
            Atemp = hess0(:,:,j);
        else
            Atemp = hessp0(:,:,j-HTHetds.nphase);
        end    
        [U11, U12, UE21, U22] = RicattiCoeffHTHet(Q0U,Atemp,HTHetds.npos);
        tmp = (U22*YU - YU*U11 + UE21 - YU*U12*YU)';   
        if j<=HTHetds.nphase
            result(HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize,HTHetds.ncoords+j) = tmp(:);
        else
            k = j-HTHetds.nphase;
            result(HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize,HTHetds.PeriodIdx+k) = tmp(:);
        end
    end
    end
    if HTHetds.YSsize    
    for j=1:HTHetds.nphase+length(HTHetds.ActiveParams)
        if j<=HTHetds.nphase
            Atemp = hess1(:,:,j);
        else
            Atemp = hessp1(:,:,j-HTHetds.nphase);
        end        
        [S11, S12, SE21, S22] = RicattiCoeffHTHet(Q0S,Atemp,HTHetds.nneg);
        tmp = (S22*YS - YS*S11 + SE21 - YS*S12*YS)';  
        result(HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+1:HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+HTHetds.YSsize,HTHetds.ncoords+HTHetds.nphase+j) = tmp(:);
    end
    end

    % UNSTABLE
    % ricatti blocks from unstable eigenspace
    if HTHetds.YUsize
    A=cjac(HTHetds.func,HTHetds.Jacobian,x0,p,HTHetds.ActiveParams);
    derU = [];
    [U11, U12, UE21, U22] = RicattiCoeffHTHet(Q0U,A,HTHetds.npos);
    U12Y = U12 * YU;
    YU12 = YU * U12;
    for j=1:HTHetds.npos
        for i=1:HTHetds.nphase-HTHetds.npos
            myres = sparse(HTHetds.nphase-HTHetds.npos,HTHetds.npos,0);
            myres(:,j) = myres(:,j) + U22(:,i);
            myres(i,:) = myres(i,:) - U11(j,:);        
            myres(i,:) = myres(i,:) - U12Y(j,:);
            myres(:,j) = myres(:,j) - YU12(:,i);        
            % Derivative is computed to YU(i,j) i=1:HTHetds.nneg j=1:HTHetds.npos
            % must come at column (j-1)*(HTHetds.nneg)+i, 
            % in form T and col-first == row-first
            rows = size(myres,1);
            cols = size(myres,2);
            for k=1:rows
                derU((k-1)*cols+1:k*cols,(j-1)*(HTHetds.nphase-HTHetds.npos)+i) = myres(k,:)';
            end
        end
    end
    end
    % STABLE
    if HTHetds.YSsize    
    A=cjac(HTHetds.func,HTHetds.Jacobian,x1,p,HTHetds.ActiveParams);
    derS = [];
    [S11, S12, SE21, S22] = RicattiCoeffHTHet(Q0S,A,HTHetds.nneg);
    S12Y = S12 * YS;
    YS12 = YS * S12;
    for i=1:HTHetds.nphase-HTHetds.nneg
        for j=1:HTHetds.nneg
            myres = sparse(HTHetds.nphase-HTHetds.nneg,HTHetds.nneg,0);
            myres(:,j) = myres(:,j) + S22(:,i);
            myres(i,:) = myres(i,:) - S11(j,:);
            myres(i,:) = myres(i,:) - S12Y(j,:);
            myres(:,j) = myres(:,j) - YS12(:,i);        
            % Derivative is computed to YS(i,j)
            % must come at column (j-1)*(HTHetds.npos)+i, 
            % in form T and col-first == row-first
            rows = size(myres,1);
            cols = size(myres,2);
            for k=1:rows
                derS((k-1)*cols+1:k*cols,(j-1)*(HTHetds.nphase-HTHetds.nneg)+i) = myres(k,:)';
            end
        end
    end
    end
    if HTHetds.YUsize    
    result(HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize,end-HTHetds.YUsize-HTHetds.YSsize+1:end-HTHetds.YSsize) = derU;
    end
    if HTHetds.YSsize    
    result(HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+1:HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+HTHetds.YSsize,end-HTHetds.YSsize+1:end) = derS;
    end

    % Component 4 (Last and first vectors along stable and unstable eigenspaces)
    % ===========
    if (HTHetds.nphase-HTHetds.npos)
        indxs = [1:HTHetds.nphase   HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase];
        vect = ups(:,1) - x0;
        Q1U = Q0U * [-YU'; eye(size(YU,1))];
        vQ = -vect' * Q0U;
        for i=1:HTHetds.nphase-HTHetds.npos
            result(HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+HTHetds.YSsize+i,indxs) = [Q1U(:,end-i+1)'     -Q1U(:,end-i+1)'];
        end
        if HTHetds.npos
            myeye = zeros(HTHetds.nphase-HTHetds.npos);
            for j=1:HTHetds.nphase-HTHetds.npos
                myeye(j,HTHetds.nphase-HTHetds.npos-j+1) = 1;
            end
            for j=1:HTHetds.npos
                result(HTHetds.ncoords+HTHetds.nphase+HTHetds.YUsize+HTHetds.YSsize+1:HTHetds.ncoords+2*HTHetds.nphase-HTHetds.npos+HTHetds.YUsize+HTHetds.YSsize,end-HTHetds.YUsize-HTHetds.YSsize+(j-1)*(HTHetds.nphase-HTHetds.npos)+1:end-HTHetds.YUsize-HTHetds.YSsize+j*(HTHetds.nphase-HTHetds.npos)) = vQ(:,j) * myeye;
            end
        end
    end

    if (HTHetds.nphase-HTHetds.nneg)
        indxs = [HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+2*HTHetds.nphase];
        vect = ups(:,end) - x1;
        Q1S = Q0S * [-YS'; eye(size(YS,1))];
        vQ = -vect' * Q0S;
        for i=1:HTHetds.nphase-HTHetds.nneg
            result(end-2-HTHetds.nphase+HTHetds.nneg+i,indxs) = [Q1S(:,end-i+1)'     -Q1S(:,end-i+1)'];
        end
        if HTHetds.nneg
            myeye = zeros(HTHetds.nphase-HTHetds.nneg);
            for j=1:HTHetds.nphase-HTHetds.nneg
                myeye(j,HTHetds.nphase-HTHetds.nneg-j+1) = 1;
            end
            for j=1:HTHetds.nneg
                result(end-2-HTHetds.nphase+HTHetds.nneg+1:end-2,end-HTHetds.YSsize+(j-1)*(HTHetds.nphase-HTHetds.nneg)+1:end-HTHetds.YSsize+j*(HTHetds.nphase-HTHetds.nneg)) = vQ(:,j) * myeye;
            end
        end
    end

    % Component 6 (Distances from endpoints to equilibrium equal to epsilons)
    % ===========
    indxs = [HTHetds.phases   HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase];
    val = ups(:,1) - x0;
    val = val' / norm(val);
    result(end-1,indxs) = [val   -val];
   
    indxs = [HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+2*HTHetds.nphase];
    val = ups(:,end) - x1;
    val = val' / norm(val);
    result(end,indxs) = [val   -val];
    result(end,end-HTHetds.YUsize-HTHetds.YSsize) = -1;

end
