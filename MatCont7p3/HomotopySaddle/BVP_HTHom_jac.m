% Build_connec2ds up jacobian of boundary value problem
%
% ============================================

function result = BVP_HTHom_jac(odefile1,x,x0,p,up,sp,T,eps1,YS,YU)
global HTHomds

if HTHomds.index == 1

    ups = reshape(x,HTHomds.nphase,HTHomds.tps);
    tmpperiod = 2 * HTHomds.T;

    % Allocate space for sparse jacobian
    result = spalloc(HTHomds.ncoords+HTHomds.npos+2,HTHomds.ncoords+HTHomds.npos+3,...
    (HTHomds.ncoords+1)^2);
    range1 = HTHomds.cols_p1;
    range2 = 1:HTHomds.ncol*HTHomds.nphase;
    range3 = 1:(HTHomds.ncol+1)*HTHomds.nphase;

    % Component 1 
    % ===========
    p = HTHomds.P0;
    p = num2cell(p);
    for i=HTHomds.tsts
        range4 = HTHomds.phases;
        % value of polynomial in each mesh point
        xval = ups(:,range1) * HTHomds.wt;
        wploc = HTHomds.wp/HTHomds.dt(i);
        % evaluate part of Jacobian matrix
        for j=HTHomds.cols
            xtmp = xval(:,j);  
        
            jac=cjac(HTHomds.func,HTHomds.Jacobian,xtmp,p,[]);
            sysjac(range4,:) = fastkron(HTHomds.ncol,HTHomds.nphase,HTHomds.wt(:,j)',jac);
        
            range4 = range4+HTHomds.nphase;
        end
        % Store result

        result(range2,range3) = [wploc-tmpperiod*sysjac];

        % shift to following intervals
        range1 = range1 + HTHomds.ncol;
        range2 = range2 + HTHomds.ncol*HTHomds.nphase;
        range3 = range3 + HTHomds.ncol*HTHomds.nphase;
    end

    % Component 2 
    % ===========
    Q0U = HTHomds.oldUnstableQ;
    if HTHomds.nneg
        for i=1:HTHomds.nneg
            result(HTHomds.ncoords-HTHomds.nphase+i,1:HTHomds.nphase) = Q0U(:,end-i+1)';    
        end
    end

    % Component 3
    % ===========
    QS = HTHomds.oldStableQ;
    if HTHomds.npos
        for i = 1:HTHomds.npos
            result(HTHomds.ncoords-HTHomds.nphase+HTHomds.nneg+i,HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords) = -1/eps1*QS(:,HTHomds.nneg+i)';  
            result(HTHomds.ncoords-HTHomds.nphase+HTHomds.nneg+i,end) = 1/eps1^2*(ups(:,end)-HTHomds.x0)'*QS(:,HTHomds.nneg+i);
        end
    end

    for i = 1:size(HTHomds.ActiveSParams,2)
        index = HTHomds.ActiveSParams(i);
        result(HTHomds.ncoords-HTHomds.nphase+HTHomds.nneg+index,HTHomds.ncoords+size(HTHomds.ActiveUParams,2)+i) = 1;
    end

    % Component 4
    % ===========
    val = ups(:,1) - HTHomds.x0;
    val = val' / norm(val);
    result(HTHomds.ncoords+1,1:HTHomds.nphase) = val;

    val = ups(:,end) - HTHomds.x0;
    val = val' / norm(val);
    result(HTHomds.ncoords+2,HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords) = val;
    
    result(HTHomds.ncoords+2,end) = -1;    

    % Component 5
    % ===========
    for k = 1:size(HTHomds.ActiveUParams,2)
        index = HTHomds.ActiveUParams(k);
        result(HTHomds.ncoords+2+index,HTHomds.ncoords+k) = HTHomds.eps0;      
    end

    for k = 1:HTHomds.npos
        result(HTHomds.ncoords+2+k,1:HTHomds.nphase) = -Q0U(:,k)';      
    end

elseif HTHomds.index == 2
   
    ups = reshape(x,HTHomds.nphase,HTHomds.tps);
    p = num2cell(p);

    pars = HTHomds.PeriodIdx + 1;
    tmpperiod = 2 * HTHomds.T;

    % Allocate space for sparse jacobian
    result = spalloc(HTHomds.ncoords+2*HTHomds.Ysize+HTHomds.nphase+2,...
        HTHomds.ncoords+HTHomds.nphase+3+2*HTHomds.Ysize,(HTHomds.ncoords+1)^2);
    range1 = HTHomds.cols_p1;
    range2 = 1:HTHomds.ncol*HTHomds.nphase;
    range3 = 1:(HTHomds.ncol+1)*HTHomds.nphase;

    % Component 1 
    % ===========
    for i=HTHomds.tsts
        range4 = HTHomds.phases;
        % value of polynomial in each mesh point
        xval = ups(:,range1) * HTHomds.wt;
        wploc = HTHomds.wp/HTHomds.dt(i);
        % evaluate part of Jacobian matrix
        for j=HTHomds.cols
            xtmp = xval(:,j);
            jac=cjac(HTHomds.func,HTHomds.Jacobian,xtmp,p,HTHomds.ActiveParams);
            sysjac(range4,:) = fastkron(HTHomds.ncol,HTHomds.nphase,HTHomds.wt(:,j)',jac);
            sysjacp(range4,:) = cjacp(HTHomds.func,HTHomds.JacobianP,xtmp,p,HTHomds.ActiveParams);
        
            range4 = range4+HTHomds.nphase;
        end
        % Store result    
        result(range2,[range3 pars]) = [wploc-tmpperiod*sysjac    -tmpperiod*sysjacp];
    
        % shift to following intervals
        range1 = range1 + HTHomds.ncol;
        range2 = range2 + HTHomds.ncol*HTHomds.nphase;
        range3 = range3 + HTHomds.ncol*HTHomds.nphase;
    end

    % Component 2 (equilibrium)
    % ===========
    result(HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords, HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.nphase) = cjac(HTHomds.func,HTHomds.Jacobian,x0,p,HTHomds.ActiveParams);
    result(HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords, HTHomds.ncoords+HTHomds.nphase+1) = cjacp(HTHomds.func,HTHomds.JacobianP,x0,p,HTHomds.ActiveParams);

    % Component 3 (Eigenspaces)
    % ===========
    Q0U = HTHomds.oldUnstableQ;
    Q0S = HTHomds.oldStableQ;
    if HTHomds.npos~=0 && HTHomds.nneg ~=0    
    hess = chess(HTHomds.func,HTHomds.Jacobian,HTHomds.Hessians,x0,p,HTHomds.ActiveParams);
    hessp = chessp(HTHomds.func,HTHomds.Jacobian,HTHomds.HessiansP,x0,p,HTHomds.ActiveParams);
    for j=1:HTHomds.nphase+1
        if j<=HTHomds.nphase
            Atemp = hess(:,:,j);
        else
            Atemp = hessp(:,:,j-HTHomds.nphase);
        end 
        [U11, U12, UE21, U22] = ricattiCoeff(Q0U,Atemp,HTHomds.npos);
        tmp = (U22*YU - YU*U11 + UE21 - YU*U12*YU)';
        result(HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.Ysize,HTHomds.ncoords+j) = tmp(:);

        [S11, S12, SE21, S22] = ricattiCoeff(Q0S,Atemp,HTHomds.nneg);
        tmp = (S22*YS - YS*S11 + SE21 - YS*S12*YS)';      
        result(HTHomds.ncoords+HTHomds.Ysize+1:HTHomds.ncoords+2*HTHomds.Ysize,HTHomds.ncoords+j) = tmp(:);    
    end

    % UNSTABLE
    % ricatti blocks from unstable eigenspace    
    A=cjac(HTHomds.func,HTHomds.Jacobian,x0,p,HTHomds.ActiveParams);
    derU = [];
    [U11, U12, UE21, U22] = ricattiCoeff(Q0U,A,HTHomds.npos);
    U12Y = U12 * YU;
    YU12 = YU * U12;
    for j=1:HTHomds.npos
        for i=1:HTHomds.nneg
            myres = sparse(HTHomds.nneg,HTHomds.npos,0);
            myres(:,j) = myres(:,j) + U22(:,i);
            myres(i,:) = myres(i,:) - U11(j,:);        
            myres(i,:) = myres(i,:) - U12Y(j,:);
            myres(:,j) = myres(:,j) - YU12(:,i);        
            % Derivative is computed to YU(i,j) i=1:HTHomds.nneg j=1:HTHomds.npos
            % must come at column (j-1)*(HTHomds.nneg)+i, 
            % in form T and col-first == row-first
            rows = size(myres,1);
            cols = size(myres,2);
            for k=1:rows
                derU((k-1)*cols+1:k*cols,(j-1)*(HTHomds.nneg)+i) = myres(k,:)';
            end
        end
    end

    % STABLE
    derS = [];
    [S11, S12, SE21, S22] = ricattiCoeff(Q0S,A,HTHomds.nneg);
    S12Y = S12 * YS;
    YS12 = YS * S12;
    for i=1:HTHomds.npos
        for j=1:HTHomds.nneg
            myres = sparse(HTHomds.npos,HTHomds.nneg,0);
            myres(:,j) = myres(:,j) + S22(:,i);
            myres(i,:) = myres(i,:) - S11(j,:);
            myres(i,:) = myres(i,:) - S12Y(j,:);
            myres(:,j) = myres(:,j) - YS12(:,i);        
            % Derivative is computed to YS(i,j)
            % must come at column (j-1)*(HTHomds.npos)+i, 
            % in form T and col-first == row-first
            rows = size(myres,1);
            cols = size(myres,2);
            for k=1:rows
                derS((k-1)*cols+1:k*cols,(j-1)*(HTHomds.npos)+i) = myres(k,:)';
            end
        end
    end
    result(HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.Ysize,end-2*HTHomds.Ysize+1:end-HTHomds.Ysize) = derU;
    result(HTHomds.ncoords+HTHomds.Ysize+1:HTHomds.ncoords+2*HTHomds.Ysize,end-HTHomds.Ysize+1:end) = derS;
    end
    
    % Component 4 (Last and first vectors along stable and unstable eigenspaces)
    % ===========
    if HTHomds.nneg
        indxs = [1:HTHomds.nphase   HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.nphase];
        vect = ups(:,1) - x0;
        Q1U = Q0U * [-YU'; eye(size(YU,1))];
        vQ = -vect' * Q0U;
        for i=1:HTHomds.nneg
            result(end-2-HTHomds.nphase+i,indxs) = [Q1U(:,end-i+1)'     -Q1U(:,end-i+1)'];
        end
        if HTHomds.npos
            myeye = zeros(HTHomds.nneg);
            for j=1:HTHomds.nneg
                myeye(j,HTHomds.nneg-j+1) = 1;
            end
            for j=1:HTHomds.npos
                result(end-2-HTHomds.nphase+1:end-2-HTHomds.nphase+HTHomds.nneg,end-2*HTHomds.Ysize+(j-1)*(HTHomds.nneg)+1:end-2*HTHomds.Ysize+j*(HTHomds.nneg)) = vQ(:,j) * myeye;
            end
        end
    end

    % Component 5 
    % ===========
    Q0S = HTHomds.oldStableQ;
    if HTHomds.npos
        Q1S = Q0S * [-YS'; eye(size(YS,1))];
        for i = 1:HTHomds.npos
            result(end-2-HTHomds.npos+i,HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords) = -1/eps1*Q1S(:,i)';    
            result(end-2-HTHomds.npos+i,HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.nphase) = 1/eps1*Q1S(:,i)';   
            result(end-2-HTHomds.npos+i,end-2*HTHomds.Ysize) = 1/eps1^2*(ups(:,end)-x0)'*Q1S(:,i);
        end
    
        if HTHomds.nneg
            vect = ups(:,end) - x0;
            vQ = vect' * Q0S;  
            myeye = zeros(HTHomds.npos);
            for j=1:HTHomds.npos
                myeye(j,j) = 1;
            end
            %afgeleide naar de elementen van YS
            for j=1:HTHomds.nneg
                result(end-2-HTHomds.npos+1:end-2,end-HTHomds.Ysize+(j-1)*HTHomds.npos+1:end-HTHomds.Ysize+j*HTHomds.npos) = vQ(:,j) * myeye;
            end
        end
    end
    result(end-2-HTHomds.npos+HTHomds.ActiveSParams,HTHomds.ncoords+HTHomds.nphase+2) = 1; %afgeleide naar actieve tau

    % Component 6 (Distances from endpoints to equilibrium equal to epsilons)
    % ===========
    indxs = [HTHomds.phases   HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.nphase];
    val = ups(:,1) - x0;
    val = val' / norm(val);
    result(end-1,indxs) = [val   -val];
    
    indxs = HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords+HTHomds.nphase;
    val = ups(:,end) - x0;
    val = val' / norm(val);
    result(end,indxs) = [val   -val];

    result(end,HTHomds.ncoords+HTHomds.nphase+3) = -1;

else 
    
    ups = reshape(x,HTHomds.nphase,HTHomds.tps);
    p = num2cell(p);

    pars = (0:1)+HTHomds.PeriodIdx +1;
    tmpperiod = 2 * T;

    % Allocate space for sparse jacobian
    result = spalloc(HTHomds.ncoords+2*HTHomds.Ysize+HTHomds.nphase+2,...
        HTHomds.ncoords+HTHomds.nphase+3+2*HTHomds.Ysize,(HTHomds.ncoords+1)^2);
    range1 = HTHomds.cols_p1;
    range2 = 1:HTHomds.ncol*HTHomds.nphase;
    range3 = 1:(HTHomds.ncol+1)*HTHomds.nphase;

    % Component 1 
    % ===========
    for i=HTHomds.tsts
        range4 = HTHomds.phases;
        % value of polynomial in each mesh point
        xval = ups(:,range1) * HTHomds.wt;
        wploc = HTHomds.wp/HTHomds.dt(i);
        % evaluate part of Jacobian matrix
        for j=HTHomds.cols
            xtmp = xval(:,j);
            frhs(range4) = feval(odefile1, 0, xtmp, p{:});
            jac=cjac(HTHomds.func,HTHomds.Jacobian,xtmp,p,HTHomds.ActiveParams);
            sysjac(range4,:) = fastkron(HTHomds.ncol,HTHomds.nphase,HTHomds.wt(:,j)',jac);
            sysjacp(range4,:) = cjacp(HTHomds.func,HTHomds.JacobianP,xtmp,p,HTHomds.ActiveParams);
        
            range4 = range4+HTHomds.nphase;
        end
        % Store result
        result(range2,[range3 pars]) = [wploc-tmpperiod*sysjac    -tmpperiod*sysjacp    -2*frhs'];

        % shift to following intervals
        range1 = range1 + HTHomds.ncol;
        range2 = range2 + HTHomds.ncol*HTHomds.nphase;
        range3 = range3 + HTHomds.ncol*HTHomds.nphase;
    end

    % Component 2 (equilibrium)
    % ===========
    result(HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords, HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.nphase) = cjac(HTHomds.func,HTHomds.Jacobian,x0,p,HTHomds.ActiveParams);
    result(HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords, HTHomds.ncoords+HTHomds.nphase+1) = cjacp(HTHomds.func,HTHomds.JacobianP,x0,p,HTHomds.ActiveParams);

    % Component 3 (Eigenspaces)
    % ===========
    Q0U = HTHomds.oldUnstableQ;
    Q0S = HTHomds.oldStableQ;
    if HTHomds.npos~=0 && HTHomds.nneg ~= 0    
    hess = chess(HTHomds.func,HTHomds.Jacobian,HTHomds.Hessians,x0,p,HTHomds.ActiveParams);
    hessp = chessp(HTHomds.func,HTHomds.Jacobian,HTHomds.HessiansP,x0,p,HTHomds.ActiveParams);
    for j=1:HTHomds.nphase+1
        if j<=HTHomds.nphase
            Atemp = hess(:,:,j);
        else
            Atemp = hessp(:,:,j-HTHomds.nphase);
        end    
        [U11, U12, UE21, U22] = ricattiCoeff(Q0U,Atemp,HTHomds.npos);
        tmp = (U22*YU - YU*U11 + UE21 - YU*U12*YU)';   
        result(HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.Ysize,HTHomds.ncoords+j) = tmp(:);
    
        [S11, S12, SE21, S22] = ricattiCoeff(Q0S,Atemp,HTHomds.nneg);
        tmp = (S22*YS - YS*S11 + SE21 - YS*S12*YS)';  
        result(HTHomds.ncoords+HTHomds.Ysize+1:HTHomds.ncoords+2*HTHomds.Ysize,HTHomds.ncoords+j) = tmp(:);
    end

    % UNSTABLE
    % ricatti blocks from unstable eigenspace
    A=cjac(HTHomds.func,HTHomds.Jacobian,x0,p,HTHomds.ActiveParams);
    derU = [];
    [U11, U12, UE21, U22] = ricattiCoeff(Q0U,A,HTHomds.npos);
    U12Y = U12 * YU;
    YU12 = YU * U12;
    for j=1:HTHomds.npos
        for i=1:HTHomds.nneg
            myres = sparse(HTHomds.nneg,HTHomds.npos,0);
            myres(:,j) = myres(:,j) + U22(:,i);
            myres(i,:) = myres(i,:) - U11(j,:);        
            myres(i,:) = myres(i,:) - U12Y(j,:);
            myres(:,j) = myres(:,j) - YU12(:,i);        
            % Derivative is computed to YU(i,j) i=1:HTHomds.nneg j=1:HTHomds.npos
            % must come at column (j-1)*(HTHomds.nneg)+i, 
            % in form T and col-first == row-first
            rows = size(myres,1);
            cols = size(myres,2);
            for k=1:rows
                derU((k-1)*cols+1:k*cols,(j-1)*(HTHomds.nneg)+i) = myres(k,:)';
            end
        end
    end

    % STABLE
    derS = [];
    [S11, S12, SE21, S22] = ricattiCoeff(Q0S,A,HTHomds.nneg);
    S12Y = S12 * YS;
    YS12 = YS * S12;
    for i=1:HTHomds.npos
        for j=1:HTHomds.nneg
            myres = sparse(HTHomds.npos,HTHomds.nneg,0);
            myres(:,j) = myres(:,j) + S22(:,i);
            myres(i,:) = myres(i,:) - S11(j,:);
            myres(i,:) = myres(i,:) - S12Y(j,:);
            myres(:,j) = myres(:,j) - YS12(:,i);        
            % Derivative is computed to YS(i,j)
            % must come at column (j-1)*(HTHomds.npos)+i, 
            % in form T and col-first == row-first
            rows = size(myres,1);
            cols = size(myres,2);
            for k=1:rows
                derS((k-1)*cols+1:k*cols,(j-1)*(HTHomds.npos)+i) = myres(k,:)';
            end
        end
    end
    result(HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.Ysize,end-2*HTHomds.Ysize+1:end-HTHomds.Ysize) = derU;
    result(HTHomds.ncoords+HTHomds.Ysize+1:HTHomds.ncoords+2*HTHomds.Ysize,end-HTHomds.Ysize+1:end) = derS;
    end

    % Component 4 (Last and first vectors along stable and unstable eigenspaces)
    % ===========
    if HTHomds.nneg
        indxs = [1:HTHomds.nphase   HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.nphase];
        vect = ups(:,1) - x0;
        Q1U = Q0U * [-YU'; eye(size(YU,1))];
        vQ = -vect' * Q0U;
        for i=1:HTHomds.nneg
            result(end-2-HTHomds.nphase+i,indxs) = [Q1U(:,end-i+1)'     -Q1U(:,end-i+1)'];
        end
        if HTHomds.npos
            myeye = zeros(HTHomds.nneg);
            for j=1:HTHomds.nneg
                myeye(j,HTHomds.nneg-j+1) = 1;
            end
            for j=1:HTHomds.npos
                result(end-2-HTHomds.nphase+1:end-2-HTHomds.nphase+HTHomds.nneg,end-2*HTHomds.Ysize+(j-1)*(HTHomds.nneg)+1:end-2*HTHomds.Ysize+j*(HTHomds.nneg)) = vQ(:,j) * myeye;
            end
        end
    end

    if HTHomds.npos
        indxs = HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords+HTHomds.nphase;
        Q1S = Q0S * [-YS'; eye(size(YS,1))];
        for i=1:HTHomds.npos
            result(end-2-HTHomds.nphase+HTHomds.nneg+i,indxs) = [Q1S(:,end-i+1)'     -Q1S(:,end-i+1)'];
        end
        if HTHomds.nneg
            vect = ups(:,end) - x0;        
            vQ = -vect' * Q0S;
            myeye = zeros(HTHomds.npos);
            for j=1:HTHomds.npos
                myeye(j,HTHomds.npos-j+1) = 1;
            end
            for j=1:HTHomds.nneg
                result(end-2-HTHomds.nphase+HTHomds.nneg+1:end-2,end-HTHomds.Ysize+(j-1)*HTHomds.npos+1:end-HTHomds.Ysize+j*HTHomds.npos) = vQ(:,j) * myeye;
            end
        end
    end

    % Component 6 (Distances from endpoints to equilibrium equal to epsilons)
    % ===========
    indxs = [HTHomds.phases   HTHomds.ncoords+1:HTHomds.ncoords+HTHomds.nphase];
    val = ups(:,1) - x0;
    val = val' / norm(val);
    result(end-1,indxs) = [val   -val];
   
    indxs = HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords+HTHomds.nphase;
    val = ups(:,end) - x0;
    val = val' / norm(val);
    result(end,indxs) = [val   -val];
    result(end,HTHomds.ncoords+HTHomds.nphase+3) = -1;

end