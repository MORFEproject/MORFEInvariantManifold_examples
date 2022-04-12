% Build_connec2ds up jacobian of boundary value problem
%
% ============================================

function result = BVP_Het_jac(odefile1,x,x0,x1,p,T,eps0,eps1,YS,YU)
global hetds

ups = reshape(x,hetds.nphase,hetds.tps);
p = num2cell(p);

np = hetds.nphase-hetds.npos-hetds.nneg+2;
if hetds.extravec(1)
    Tvar = 1;
    pars = (0:np)+hetds.PeriodIdx +1;
else
    Tvar = 0;
    pars = (0:(np-1))+hetds.PeriodIdx +1;
end
tmpperiod = 2 * T;

phasepresent = 0;
if sum(hetds.extravec) == 2
    phasepresent = 1;
elseif sum(hetds.extravec) ~= 1
    error('Wrong number of free heteroclinic parameters');
end
    
% Allocate space for sparse jacobian
result = spalloc(hetds.ncoords+hetds.YUsize+hetds.YSsize+3*hetds.nphase-hetds.nneg-hetds.npos+phasepresent+2,...
hetds.ncoords+3*hetds.nphase-hetds.npos-hetds.nneg+sum(hetds.extravec)+2+hetds.YUsize+hetds.YSsize,(hetds.ncoords+1)^2);
range1 = hetds.cols_p1;
range2 = 1:hetds.ncol*hetds.nphase;
range3 = 1:(hetds.ncol+1)*hetds.nphase;

% Component 1 
% ===========
for i=hetds.tsts
    range4 = hetds.phases;
    % value of polynomial in each mesh point
    xval = ups(:,range1) * hetds.wt;
    wploc = hetds.wp/hetds.dt(i);
    % evaluate part of Jacobian matrix
    for j=hetds.cols
        xtmp = xval(:,j);
        if Tvar
            frhs(range4) = feval(odefile1, 0, xtmp, p{:});
        end
        jac=cjac(hetds.func,hetds.Jacobian,xtmp,p,hetds.ActiveParams);
        sysjac(range4,:) = fastkron(hetds.ncol,hetds.nphase,hetds.wt(:,j)',jac);
        sysjacp(range4,:) = cjacp(hetds.func,hetds.JacobianP,xtmp,p,hetds.ActiveParams);
        range4 = range4+hetds.nphase;
    end
    % Store result
    if Tvar
        result(range2,[range3 pars]) = [wploc-tmpperiod*sysjac    -tmpperiod*sysjacp    -2*frhs'];
    else
        result(range2,[range3 pars]) = [wploc-tmpperiod*sysjac    -tmpperiod*sysjacp];
    end
    % shift to following intervals
    range1 = range1 + hetds.ncol;
    range2 = range2 + hetds.ncol*hetds.nphase;
    range3 = range3 + hetds.ncol*hetds.nphase;
end

% Component 2 (equilibrium)
% ===========
result(hetds.ncoords-hetds.nphase+1:hetds.ncoords, hetds.ncoords+1:hetds.ncoords+hetds.nphase) = cjac(hetds.func,hetds.Jacobian,x0,p,hetds.ActiveParams);
result(hetds.ncoords-hetds.nphase+1:hetds.ncoords, hetds.PeriodIdx+1:hetds.PeriodIdx+length(hetds.ActiveParams)) = cjacp(hetds.func,hetds.JacobianP,x0,p,hetds.ActiveParams);
result(hetds.ncoords+1:hetds.ncoords+hetds.nphase, hetds.ncoords+hetds.nphase+1:hetds.ncoords+2*hetds.nphase) = cjac(hetds.func,hetds.Jacobian,x1,p,hetds.ActiveParams);
result(hetds.ncoords+1:hetds.ncoords+hetds.nphase, hetds.PeriodIdx+1:hetds.PeriodIdx+length(hetds.ActiveParams)) = cjacp(hetds.func,hetds.JacobianP,x1,p,hetds.ActiveParams);    

% Component 3 (phase condition)
% ===========
if phasepresent
    range1 = hetds.cols_p1;
    range3 = 1:(hetds.ncol+1)*hetds.nphase;
    icjac = zeros(1,hetds.ncoords);
    for i=1:hetds.ntst
        tmp = hetds.dt(i)*(hetds.upoldp(:,range1) .* hetds.pwi);
        icjac(range3) = icjac(range3) + tmp(1:(hetds.ncol+1)*hetds.nphase);
        range1 = range1 + hetds.ncol;
        range3 = range3 + hetds.ncol*hetds.nphase;
    end
    result(hetds.ncoords+1,[hetds.coords]) = icjac;
end

% Component 4 (Eigenspaces)
% ===========
Q0U = hetds.oldUnstableQ;
Q0S = hetds.oldStableQ;
hess0 = chess(hetds.func,hetds.Jacobian,hetds.Hessians,x0,p,hetds.ActiveParams);
hessp0 = chessp(hetds.func,hetds.Jacobian,hetds.HessiansP,x0,p,hetds.ActiveParams);
hess1 = chess(hetds.func,hetds.Jacobian,hetds.Hessians,x1,p,hetds.ActiveParams);
hessp1 = chessp(hetds.func,hetds.Jacobian,hetds.HessiansP,x1,p,hetds.ActiveParams);
for j=1:hetds.nphase+length(hetds.ActiveParams)
    if j<=hetds.nphase
        Atemp = hess0(:,:,j);
    else
        Atemp = hessp0(:,:,j-hetds.nphase);
    end
    [U11, U12, UE21, U22] = RicattiCoeffHet(Q0U,Atemp,hetds.npos);
    tmp = (U22*YU - YU*U11 + UE21 - YU*U12*YU)';   
    if j<=hetds.nphase
        result(hetds.ncoords+hetds.nphase+phasepresent+1:hetds.ncoords+hetds.nphase+phasepresent+hetds.YUsize,hetds.ncoords+j) = tmp(:);
    else
        k = j-hetds.nphase;  
        result(hetds.ncoords+hetds.nphase+phasepresent+1:hetds.ncoords+hetds.nphase+phasepresent+hetds.YUsize,hetds.PeriodIdx+k) = tmp(:); 
    end
end

for j=1:hetds.nphase+length(hetds.ActiveParams)
    if j<=hetds.nphase
        Atemp = hess1(:,:,j);
    else
        Atemp = hessp1(:,:,j-hetds.nphase);
    end
    [S11, S12, SE21, S22] = RicattiCoeffHet(Q0S,Atemp,hetds.nneg);
    tmp = (S22*YS - YS*S11 + SE21 - YS*S12*YS)';  
    result(hetds.ncoords+hetds.nphase+phasepresent+hetds.YUsize+1:hetds.ncoords+hetds.nphase+phasepresent+hetds.YUsize+hetds.YSsize,hetds.ncoords+hetds.nphase+j) = tmp(:);
end

% UNSTABLE
% ricatti blocks from unstable eigenspace
A=cjac(hetds.func,hetds.Jacobian,x0,p,hetds.ActiveParams);
derU = [];
[U11, U12, UE21, U22] = RicattiCoeffHet(Q0U,A,hetds.npos);
U12Y = U12 * YU;
YU12 = YU * U12;
for j=1:hetds.npos
    for i=1:hetds.nphase-hetds.npos
        myres = sparse(hetds.nphase-hetds.npos,hetds.npos,0);
        myres(:,j) = myres(:,j) + U22(:,i);
        myres(i,:) = myres(i,:) - U11(j,:);        
        myres(i,:) = myres(i,:) - U12Y(j,:);
        myres(:,j) = myres(:,j) - YU12(:,i);        
        % Derivative is computed to YU(i,j) i=1:hetds.nneg j=1:hetds.npos
        % must come at column (j-1)*(hetds.nneg)+i, 
        % in form T and col-first == row-first
        rows = size(myres,1);
        cols = size(myres,2);   
        for k=1:rows    
            derU((k-1)*cols+1:k*cols,(j-1)*(hetds.nphase-hetds.npos)+i) = myres(k,:)';
        end
    end
end

% STABLE
A=cjac(hetds.func,hetds.Jacobian,x1,p,hetds.ActiveParams);
derS = [];
[S11, S12, SE21, S22] = RicattiCoeffHet(Q0S,A,hetds.nneg);
S12Y = S12 * YS;
YS12 = YS * S12;
for i=1:hetds.nphase-hetds.nneg
    for j=1:hetds.nneg
        myres = sparse(hetds.nphase-hetds.nneg,hetds.nneg,0);
        myres(:,j) = myres(:,j) + S22(:,i);
        myres(i,:) = myres(i,:) - S11(j,:);
        myres(i,:) = myres(i,:) - S12Y(j,:); 
        myres(:,j) = myres(:,j) - YS12(:,i);        
        % Derivative is computed to YS(i,j)
        % must come at column (j-1)*(hetds.npos)+i, 
        % in form T and col-first == row-first
        rows = size(myres,1);
        cols = size(myres,2);
        for k=1:rows
            derS((k-1)*cols+1:k*cols,(j-1)*(hetds.nphase-hetds.nneg)+i) = myres(k,:)';
        end
    end
end
result(hetds.ncoords+hetds.nphase+phasepresent+1:hetds.ncoords+hetds.nphase+phasepresent+hetds.YUsize,end-hetds.YUsize-hetds.YSsize+1:end-hetds.YSsize) = derU;
result(hetds.ncoords+hetds.nphase+phasepresent+hetds.YUsize+1:hetds.ncoords+hetds.nphase+phasepresent+hetds.YUsize+hetds.YSsize,end-hetds.YSsize+1:end) = derS;

% Component 5 (Last and first vectors along stable and unstable
% eigenspaces)
% ===========
if (hetds.nphase-hetds.npos)
    indxs = [1:hetds.nphase   hetds.ncoords+1:hetds.ncoords+hetds.nphase];
    vect = ups(:,1) - x0;
    Q1U = Q0U * [-YU'; eye(size(YU,1))];
    vQ = -vect' * Q0U;    
    for i=1:hetds.nphase-hetds.npos
        result(hetds.ncoords+hetds.nphase+phasepresent+hetds.YUsize+hetds.YSsize+i,indxs) = [Q1U(:,end-i+1)'     -Q1U(:,end-i+1)'];
    end
    myeye = zeros(hetds.nphase-hetds.npos);    
    for j=1:hetds.nphase-hetds.npos    
        myeye(j,hetds.nphase-hetds.npos-j+1) = 1;        
    end
    for j=1:hetds.npos    
        result(hetds.ncoords+hetds.nphase+phasepresent+hetds.YUsize+hetds.YSsize+1:hetds.ncoords+phasepresent+2*hetds.nphase-hetds.npos+hetds.YUsize+hetds.YSsize,end-hetds.YUsize-hetds.YSsize+(j-1)*(hetds.nphase-hetds.npos)+1:end-hetds.YUsize-hetds.YSsize+j*(hetds.nphase-hetds.npos)) = vQ(:,j) * myeye;                
    end
end

if (hetds.nphase-hetds.nneg)
    indxs = [hetds.ncoords-hetds.nphase+1:hetds.ncoords hetds.ncoords+hetds.nphase+1:hetds.ncoords+2*hetds.nphase];
    vect = ups(:,end) - x1;
    Q1S = Q0S * [-YS'; eye(size(YS,1))];
    vQ = -vect' * Q0S;   
    for i=1:(hetds.nphase-hetds.nneg)         
        result(end-2-hetds.nphase+hetds.nneg+i,indxs) = [Q1S(:,end-i+1)'     -Q1S(:,end-i+1)'];        
    end
    myeye = zeros(hetds.nphase-hetds.nneg);    
    for j=1:hetds.nphase-hetds.nneg        
        myeye(j,hetds.nphase-hetds.nneg-j+1) = 1;                
    end
    for j=1:hetds.nneg    
        result(end-2-hetds.nphase+hetds.nneg+1:end-2,end-hetds.YSsize+(j-1)*(hetds.nphase-hetds.nneg)+1:end-hetds.YSsize+j*(hetds.nphase-hetds.nneg)) = vQ(:,j) * myeye;               
    end
end

% Component 6 (Distances from endpoints to equilibrium equal to epsilons)
% ===========
indxs = [hetds.phases   hetds.ncoords+1:hetds.ncoords+hetds.nphase];
val = ups(:,1) - x0;
val = val' / norm(val);
result(end-1,indxs) = [val   -val];
if hetds.extravec(2)
    result(end-1,end-hetds.YUsize-hetds.YSsize-hetds.extravec(3)) = -1;
end

indxs = [hetds.ncoords-hetds.nphase+1:hetds.ncoords hetds.ncoords+hetds.nphase+1:hetds.ncoords+2*hetds.nphase];
val = ups(:,end) - x1;
val = val' / norm(val);
result(end,indxs) = [val   -val];
if hetds.extravec(3)
    result(end,end-hetds.YUsize-hetds.YSsize) = -1;
end
