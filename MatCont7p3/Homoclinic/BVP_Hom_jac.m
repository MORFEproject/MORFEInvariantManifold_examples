% Buihomds up jacobian of boundary value problem
%
% ============================================

function result = BVP_Hom_jac(odefile1,x,x0,p,T,eps0,eps1,YS,YU)
global homds xvar
ups = reshape(x,homds.nphase,homds.tps);
p = num2cell(p);

xvar = x;

if homds.extravec(1)
    Tvar = 1;
    pars = (0:2)+homds.PeriodIdx +1;
else
    Tvar = 0;
    pars = (0:1)+homds.PeriodIdx +1;
end
tmpperiod = 2 * T;

phasepresent = 0;
if sum(homds.extravec) == 2
    phasepresent = 1;
elseif sum(homds.extravec) ~= 1
    error('Wrong number of free homoclinic parameters');
end

% Allocate space for sparse jacobian
result = spalloc(homds.ncoords+phasepresent+2*homds.Ysize+homds.nphase+2,...
    homds.ncoords+homds.nphase+2+sum(homds.extravec)+2*homds.Ysize,(homds.ncoords+1)^2);
range1 = homds.cols_p1;
range2 = 1:homds.ncol*homds.nphase;
range3 = 1:(homds.ncol+1)*homds.nphase;

% Component 1 
% ===========
for i=homds.tsts
    range4 = homds.phases;
    % value of polynomial in each mesh point
    xval = ups(:,range1) * homds.wt;
    wploc = homds.wp/homds.dt(i);
    % evaluate part of Jacobian matrix
    for j=homds.cols
        xtmp = xval(:,j);
        if Tvar
            frhs(range4) = feval(odefile1, 0, xtmp, p{:});
        end
        jac=cjac(homds.func,homds.Jacobian,xtmp,p,homds.ActiveParams);
        sysjac(range4,:) = fastkron(homds.ncol,homds.nphase,homds.wt(:,j)',jac);
        sysjacp(range4,:) = cjacp(homds.func,homds.JacobianP,xtmp,p,homds.ActiveParams);
        
        range4 = range4+homds.nphase;
    end
    % Store result
    if Tvar && (homds.HTPstep == 0)
        result(range2,[range3 pars]) = [wploc-tmpperiod*sysjac    -tmpperiod*sysjacp    -2*frhs'];
    else
        result(range2,[range3 pars]) = [wploc-tmpperiod*sysjac    -tmpperiod*sysjacp];
    end    
    % shift to following intervals
    range1 = range1 + homds.ncol;
    range2 = range2 + homds.ncol*homds.nphase;
    range3 = range3 + homds.ncol*homds.nphase;
end

% Component 2 (equilibrium)
% ===========
result(homds.ncoords-homds.nphase+1:homds.ncoords, homds.ncoords+1:homds.ncoords+homds.nphase) = cjac(homds.func,homds.Jacobian,x0,p,homds.ActiveParams);
result(homds.ncoords-homds.nphase+1:homds.ncoords, homds.ncoords+homds.nphase+1:homds.ncoords+homds.nphase+2) = cjacp(homds.func,homds.JacobianP,x0,p,homds.ActiveParams);

% Component 3 (phase condition)
% ===========
if phasepresent
    range1 = homds.cols_p1;
    range3 = 1:(homds.ncol+1)*homds.nphase;
    icjac = zeros(1,homds.ncoords);
    for i=1:homds.ntst
        tmp = homds.dt(i)*(homds.upoldp(:,range1) .* homds.pwi);
        icjac(range3) = icjac(range3) + tmp(1:(homds.ncol+1)*homds.nphase);
        range1 = range1 + homds.ncol;
        range3 = range3 + homds.ncol*homds.nphase;
    end
    result(homds.ncoords+1,[homds.coords]) = icjac;
end

% Component 4 (Eigenspaces)
% ===========
if homds.HTPstep == 0
Q0U = homds.oldUnstableQ;
Q0S = homds.oldStableQ;
hess = chess(homds.func,homds.Jacobian,homds.Hessians,x0,p,homds.ActiveParams);
hessp = chessp(homds.func,homds.Jacobian,homds.HessiansP,x0,p,homds.ActiveParams);
for j=1:homds.nphase+2
    if j<=homds.nphase
        Atemp = hess(:,:,j);
    else
        Atemp = hessp(:,:,j-homds.nphase);
    end    
    [U11, U12, UE21, U22] = ricattiCoeff(Q0U,Atemp,homds.npos);
    tmp = (U22*YU - YU*U11 + UE21 - YU*U12*YU)';
   
  %  for i=1:homds.nneg
        result(homds.ncoords+phasepresent+1:homds.ncoords+phasepresent+homds.Ysize,homds.ncoords+j) = tmp(:);
  % end    
    [S11, S12, SE21, S22] = ricattiCoeff(Q0S,Atemp,homds.nneg);
    tmp = (S22*YS - YS*S11 + SE21 - YS*S12*YS)';
  %  for i=1:homds.npos
        result(homds.ncoords+phasepresent+homds.Ysize+1:homds.ncoords+phasepresent+2*homds.Ysize,homds.ncoords+j) = tmp(:);
  %  end
end

% UNSTABLE
% ricatti blocks from unstable eigenspace
A=cjac(homds.func,homds.Jacobian,x0,p,homds.ActiveParams);
derU = [];
[U11, U12, UE21, U22] = ricattiCoeff(Q0U,A,homds.npos);
U12Y = U12 * YU;
YU12 = YU * U12;
for j=1:homds.npos
    for i=1:homds.nneg
        myres = sparse(homds.nneg,homds.npos,0);
        myres(:,j) = myres(:,j) + U22(:,i);
        myres(i,:) = myres(i,:) - U11(j,:);        
        myres(i,:) = myres(i,:) - U12Y(j,:);
        myres(:,j) = myres(:,j) - YU12(:,i);        
        % Derivative is computed to YU(i,j) i=1:homds.nneg j=1:homds.npos
        % must come at column (j-1)*(homds.nneg)+i, 
        % in form T and col-first == row-first
        rows = size(myres,1);
        cols = size(myres,2);
        for k=1:rows
            derU((k-1)*cols+1:k*cols,(j-1)*(homds.nneg)+i) = myres(k,:)';
        end
    end
end

% STABLE
derS = [];
[S11, S12, SE21, S22] = ricattiCoeff(Q0S,A,homds.nneg);
S12Y = S12 * YS;
YS12 = YS * S12;
for i=1:homds.npos
    for j=1:homds.nneg
        myres = sparse(homds.npos,homds.nneg,0);
        myres(:,j) = myres(:,j) + S22(:,i);
        myres(i,:) = myres(i,:) - S11(j,:);
        myres(i,:) = myres(i,:) - S12Y(j,:);
        myres(:,j) = myres(:,j) - YS12(:,i);        
        % Derivative is computed to YS(i,j)
        % must come at column (j-1)*(homds.npos)+i, 
        % in form T and col-first == row-first
        rows = size(myres,1);
        cols = size(myres,2);
        for k=1:rows
            derS((k-1)*cols+1:k*cols,(j-1)*(homds.npos)+i) = myres(k,:)';
        end
    end
end
result(homds.ncoords+phasepresent+1:homds.ncoords+phasepresent+homds.Ysize,end-2*homds.Ysize+1:end-homds.Ysize) = derU;
result(homds.ncoords+phasepresent+homds.Ysize+1:homds.ncoords+phasepresent+2*homds.Ysize,end-homds.Ysize+1:end) = derS;
end

% Component 5 (Last and first vectors along stable and unstable eigenspaces)
% ===========
if homds.npos
    indxs = [1:homds.nphase   homds.ncoords+1:homds.ncoords+homds.nphase];
    vect = ups(:,1) - x0;
    Q1U = Q0U * [-YU'; eye(size(YU,1))];
    vQ = -vect' * Q0U;
    for i=1:homds.nneg
        result(end-2-homds.nphase+i,indxs) = [Q1U(:,end-i+1)'     -Q1U(:,end-i+1)'];
    end
    myeye = zeros(homds.nneg);
    for j=1:homds.nneg
        myeye(j,homds.nneg-j+1) = 1;
    end
    for j=1:homds.npos
        result(end-2-homds.nphase+1:end-2-homds.nphase+homds.nneg,end-2*homds.Ysize+(j-1)*(homds.nneg)+1:end-2*homds.Ysize+j*(homds.nneg)) = vQ(:,j) * myeye;
    end
end

if homds.nneg
    indxs = homds.ncoords-homds.nphase+1:homds.ncoords+homds.nphase;
    vect = ups(:,end) - x0;
    Q1S = Q0S * [-YS'; eye(size(YS,1))];
    vQ = -vect' * Q0S;
    for i=1:homds.npos
        result(end-2-homds.nphase+homds.nneg+i,indxs) = [Q1S(:,end-i+1)'     -Q1S(:,end-i+1)'];
    end
    myeye = zeros(homds.npos);
    for j=1:homds.npos
        myeye(j,homds.npos-j+1) = 1;
    end
    for j=1:homds.nneg
        result(end-2-homds.nphase+homds.nneg+1:end-2,end-homds.Ysize+(j-1)*homds.npos+1:end-homds.Ysize+j*homds.npos) = vQ(:,j) * myeye;
    end
end

% Component 6 (Distances from endpoints to equilibrium equal to epsilons)
% ===========
    indxs = [homds.phases   homds.ncoords+1:homds.ncoords+homds.nphase];
    val = ups(:,1) - x0;
    val = val' / norm(val);
    result(end-1,indxs) = [val   -val];
   
    if homds.extravec(2)
        result(end-1,homds.ncoords+homds.nphase+3+homds.extravec(1)) = -1;
    end
    
    indxs = homds.ncoords-homds.nphase+1:homds.ncoords+homds.nphase;
    val = ups(:,end) - x0;
    val = val' / norm(val);
    result(end,indxs) = [val   -val];
    if homds.extravec(3)
        result(end,homds.ncoords+homds.nphase+3+sum(homds.extravec(1:2))) = -1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
