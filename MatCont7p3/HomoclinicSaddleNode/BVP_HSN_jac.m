% Builds up jacobian of boundary value problem
%
% ============================================

function result = BVP_HSN_jac(odefile1,x,x0,p,T,eps0,eps1,YS,YU)

global homds

ups = reshape(x,homds.nphase,homds.tps);
p = num2cell(p);

if homds.extravec(1)
    Tvar = 1;
    pars = (0:2)+homds.PeriodIdx +1;
else
    Tvar = 0;
    pars = (0:1)+homds.PeriodIdx +1;
end
tmpperiod = 2 * T;

% Allocate space for sparse jacobian
result = spalloc(homds.ncoords+1+homds.Ysize+homds.nphase+2,...
    homds.ncoords+homds.nphase+2+sum(homds.extravec)+homds.Ysize,(homds.ncoords+1)^2);
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
        jac = cjac(homds.func,homds.Jacobian,xtmp,p,homds.ActiveParams);
        sysjac(range4,:) = fastkron(homds.ncol,homds.nphase,homds.wt(:,j)',jac);
        sysjacp(range4,:) = cjacp(homds.func,homds.JacobianP,xtmp,p,homds.ActiveParams);
        
        range4 = range4+homds.nphase;
    end
    % Store result
    if Tvar 
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
tmpjac = cjac(homds.func,homds.Jacobian,x0,p,homds.ActiveParams);
result(homds.ncoords-homds.nphase+1:homds.ncoords, homds.ncoords+1:homds.ncoords+homds.nphase) = tmpjac;
tmp1 = cjacp(homds.func,homds.JacobianP,x0,p,homds.ActiveParams);
result(homds.ncoords-homds.nphase+1:homds.ncoords,  homds.ncoords+homds.nphase+1:homds.ncoords+homds.nphase+2) = tmp1;

% Component 2bis (limit point)
% ==============
  Bord = [tmpjac homds.wvector;...
        homds.vvector' 0];
  bunit = [zeros(homds.nphase,1);1];
  vext = Bord \ bunit;
  wext = Bord' \ bunit;
  jac = [jac tmp1]; 
  hess = chess(homds.func,homds.Jacobian,homds.Hessians,x0,p,homds.ActiveParams);
  for i=1:homds.nphase
      jac(homds.nphase+1,i)=-wext(1:homds.nphase)'*hess(:,:,i)*vext(1:homds.nphase);
  end
  hessp = chessp(homds.func,homds.Jacobian,homds.HessiansP,x0,p,homds.ActiveParams);
  for i=1:length(homds.ActiveParams);
      jac(homds.nphase+1,homds.nphase+i)=-wext(1:homds.nphase)'*hessp(:,:,i)*vext(1:homds.nphase);
  end
  result(homds.ncoords+1,homds.ncoords+1:homds.ncoords+homds.nphase+2) = jac(end,:);

% Component 3 (phase condition)
% ===========
range1 = homds.cols_p1;
range3 = 1:(homds.ncol+1)*homds.nphase;
icjac = zeros(1,homds.ncoords);
for i=1:homds.ntst
    tmp = homds.dt(i)*(homds.upoldp(:,range1) .* homds.pwi);
    icjac(range3) = icjac(range3) + tmp(1:(homds.ncol+1)*homds.nphase);
    
    range1 = range1 + homds.ncol;
    range3 = range3 + homds.ncol*homds.nphase;
end
result(homds.ncoords+2,[homds.coords]) = icjac;


% Component 4 (Eigenspaces)
% ===========
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
    if ~isempty(YU)
        [U11, U12, UE21, U22] = ricattiCoeff(Q0U,Atemp,homds.npos);
        tmp = (U22*YU - YU*U11 + UE21 - YU*U12*YU)';
        result(homds.ncoords+2+1:homds.ncoords+2+homds.npos*(homds.nneg+1),homds.ncoords+j) = tmp(:);
    end
    if ~isempty(YS)
        [S11, S12, SE21, S22] = ricattiCoeff(Q0S,Atemp,homds.nneg);
        tmp = (S22*YS - YS*S11 + SE21 - YS*S12*YS)';
        result(homds.ncoords+2+homds.npos*(homds.nneg+1)+1:homds.ncoords+2+homds.npos*(homds.nneg+1)+homds.nneg*(homds.npos+1),homds.ncoords+j) = tmp(:);
    end
end

% UNSTABLE
% Riccati blocks from unstable eigenspace
derU = [];
if ~isempty(YU)
    [U11, U12, UE21, U22] = ricattiCoeff(Q0U,tmpjac,homds.npos);
    U12Y = U12 * YU;
    YU12 = YU * U12;
end
for j=1:homds.npos
    for i=1:homds.nneg+1
        myres = sparse(homds.nneg+1,homds.npos,0);
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
            derU((k-1)*cols+1:k*cols,(j-1)*(homds.nneg+1)+i) = myres(k,:)';
        end
    end
end

% STABLE
derS = [];
if ~isempty(YS)
    [S11, S12, SE21, S22] = ricattiCoeff(Q0S,tmpjac,homds.nneg);
    S12Y = S12 * YS;
    YS12 = YS * S12;
end
for i=1:homds.npos+1
    for j=1:homds.nneg
        myres = sparse(homds.npos+1,homds.nneg,0);
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
            derS((k-1)*cols+1:k*cols,(j-1)*(homds.npos+1)+i) = myres(k,:)';
        end
    end
end
if ~isempty(YU)
    result(homds.ncoords+2+1:homds.ncoords+2+homds.npos*(homds.nneg+1),end-homds.Ysize+1:end-(homds.npos+1)*(homds.nneg)) = derU;
end
if ~isempty(YS)
    result(homds.ncoords+2+homds.npos*(homds.nneg+1)+1:homds.ncoords+2+homds.Ysize,end-(homds.npos+1)*(homds.nneg)+1:end) = derS;
end
    

% Component 5 (Last and first vectors along stable and unstable eigenspaces)
% ===========
opt.disp = 0;
[V,D] = eig(tmpjac);
[Y,i] = min(abs(diag(D)));
V = V(:,i);
[W,D] = eig(tmpjac');
[Y,i] = min(abs(diag(D)));
W = W(:,i);
D = D(i,i);

almostA = tmpjac - D * eye(size(tmpjac));
bigMat = [almostA W; V' 0];
hessx = chess(homds.func,homds.Jacobian,homds.Hessians,x0,p,homds.ActiveParams);
hessp = chessp(homds.func,homds.Jacobian,homds.HessiansP,x0,p,homds.ActiveParams);
for i=1:homds.nphase
    dlambdadx = (W' * hessx(:,:,i) * V) / (W' * V);
    rhs = -hessx(:,:,i) * V + dlambdadx * V;
    dVdxg = bigMat \ [rhs ;0];
    dV(:,i) = dVdxg(1:end-1,1);
end
for i=1:length(homds.ActiveParams)
    dlambdadp = (W' * hessp(:,:,i)  * V) / (W' * V);
    rhs2 = -hessp(:,:,i) * V + dlambdadp * V;
    dVdpg = bigMat \ [rhs2 ;0];
    dV(:,homds.nphase+i) = dVdpg(1:end-1,1);
end

if homds.nneg
    indxs = [1:homds.nphase   homds.ncoords+1:homds.ncoords+homds.nphase];
    vect = ups(:,1) - x0;
    if homds.npos
        Q1U = Q0U * [eye(size(YU,1)); YU];
    else
        Q1U = [];
    end;
    
    [Q,R] = qr(V);
    for j=1:homds.nneg
        dQtot{j} = [];
    end
    for i=1:homds.nphase+2
        Qdv = (Q' * dV(:,i));
        Qdv = tril(Qdv,-1);
        B = Qdv / R;
        B = tril(B) - tril(B)';
        dQ = Q * B;
        for j=1:homds.nneg            
            tmp = dQtot{j};
            dQtot{j} = [tmp dQ(:,end-j+1)];
        end
    end    
    
     myA = [Q1U V];
     [Qn,Rn] = qr(myA);
     
    vQ = -vect' * [Q0U V];
    
    for i=1:homds.nneg
        result(end-1-homds.nphase+i,indxs) = [Qn(:,end-i+1)'     -Qn(:,end-i+1)'];
    end
    for i=1:homds.nneg
        vQb = vect' * dQtot{i};
        result(end-1-homds.nphase+i,homds.ncoords+1:homds.ncoords+homds.nphase+2) = ...
            result(end-1-homds.nphase+i,homds.ncoords+1:homds.ncoords+homds.nphase+2) + ...
            vQb;
    end
    if homds.npos
        for j=1:homds.npos
            result(end-homds.nphase:end-homds.nphase-1+homds.nneg,end-homds.Ysize+(j-1)*(homds.nneg+1)+1:end-homds.Ysize+j*(homds.nneg+1)) = vQ(:,j) .* eye(homds.nneg,homds.nneg+1);
        end
    end
end

if homds.npos
    indxs = [homds.ncoords-homds.nphase+1:homds.ncoords+homds.nphase];
    vect = ups(:,end) - x0;
    if homds.nneg
        Q1S = Q0S * [eye(size(YS,1)); YS];
    else 
        Q1S = [];
    end
    [Q,R] = qr(V);
    for j=1:homds.npos
        dQtot{j} = [];
    end
    for i=1:homds.nphase+2
        Qdv = (Q' * dV(:,i));
        Qdv = tril(Qdv,-1);
        B = Qdv / R;
        B = tril(B) - tril(B)';
        dQ = Q * B;
        for j=1:homds.npos            
            tmp = dQtot{j};
            dQtot{j} = [tmp dQ(:,end-j+1)];
        end
    end   
    
    myA = [Q1S V];
    [Qn,Rn] = qr(myA);
    vQ = -vect' * Q0S;    
    
    for i=1:homds.npos
        result(end-2-homds.npos+i,indxs) = [Qn(:,end-i+1)'     -Qn(:,end-i+1)'];
    end
    for i=1:homds.npos
        vQb = vect' * dQtot{i};
        result(end-2-homds.npos+i,homds.ncoords+1:homds.ncoords+homds.nphase+2) = ...
            result(end-2-homds.npos+i,homds.ncoords+1:homds.ncoords+homds.nphase+2) + ...
           vQb;
    end    
    if homds.nneg
        for j=1:homds.nneg
            result(end-2-homds.npos+1:end-2,end-(homds.npos+1)*(homds.nneg)+(j-1)*(homds.npos+1)+1:end-(homds.npos+1)*(homds.nneg)+j*(homds.npos+1)) = vQ(:,j) * eye(homds.npos,homds.npos+1);
        end
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
    
    indxs = [homds.ncoords-homds.nphase+1:homds.ncoords+homds.nphase];
    val = ups(:,end) - x0;
    val = val' / norm(val);
    result(end,indxs) = [val   -val];
    if homds.extravec(3)
        result(end,homds.ncoords+homds.nphase+3+sum(homds.extravec(1:2))) = -1;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

