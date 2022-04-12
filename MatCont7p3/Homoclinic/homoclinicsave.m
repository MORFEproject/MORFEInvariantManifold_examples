function out = homoclinic
%
% homoclinic curve definition file for a problem in odefile
% 
global homds cds
    out{1}  = @curve_func;
    out{2}  = @defaultprocessor;
    out{3}  = @options;
    out{4}  = @jacobian;
    out{5}  = @hessians;
    out{6}  = @testf;
    out{7}  = @userf;
    out{8}  = @process;
    out{9}  = @singmat;
    out{10} = @locate;
    out{11} = @init;
    out{12} = @done;
    out{13} = @adapt;
return


%----------------------------------------------------

function func = curve_func(arg)

  [x,x0,p,T,eps0,eps1,YS,YU] = rearr(arg);
  func = BVP_Hom(x,x0,p,T,eps0,eps1,YS,YU);
  
%------------------------------------------------------

function varargout = jacobian(varargin)
  global homds BigJac
  [x,x0,p,T,eps0,eps1,YS,YU] = rearr(varargin{1});
  varargout{1} = BVP_Hom_jac(homds.func,x,x0,p,T,eps0,eps1,YS,YU);
  BigJac = varargout{1};
%-----------------------------------------------------

function hessians(varargin)

%------------------------------------------------------

function varargout = defaultprocessor(varargin)
global homds cds
    
  [x,x0,p,T,eps0,eps1,YS,YU] = rearr(varargin{1});
  v = rearr(varargin{2});
  
  homds.ndim = length(varargin{1});

  % update
  if ~isempty(homds.ups)
    homds.upold = homds.ups;
  end
  homds.ups = reshape(x,homds.nphase,homds.tps);
  homds.vps = reshape(v,homds.nphase,homds.tps);
  
  % update upoldp
  p1 = num2cell(p);
  for i=1:homds.tps
    homds.upoldp(:,i) = 2*T*feval(homds.func, 0, homds.ups(:,i), p1{:});
  end
  homds.eps0 = eps0;
  homds.eps1 = eps1;
  homds.YS = YS;
  homds.YU = YU;
  homds.T = T;
  homds.x0 = x0;
  
% Update dimensions
% -----------------
p = num2cell(p);
A = cjac(homds.func,homds.Jacobian,x0,p,homds.ActiveParams);
D = eig(A);
% nneg = dimension of stable subspace
homds.nneg = sum(real(D) < 0);

% If one eigenvalue is (practically) zero, and the one of the subspaces has
% zero dimension, change this dimension with 1.
if (homds.nneg == homds.nphase)
    if min(abs(real(D))) < 1e-2
        homds.nneg = homds.nneg -1;
    end
end
if (homds.nneg == 0)
    if min(abs(real(D))) < 1e-2
        homds.nneg = homds.nneg +1;
    end
end
homds.npos = homds.nphase-homds.nneg;
homds.Ysize = homds.nneg*homds.npos;
  
  if nargin > 2
    % set data in special point structure
    s = varargin{3};
    s.data.timemesh = homds.msh;
    s.data.ntst = homds.ntst;
    s.data.ncol = homds.ncol;
    s.data.parametervalues = p;
    s.data.T = T;
    varargout{3} = s;
  end
  % all done succesfully
  varargout{1} = 0;
  varargout{2} = homds.msh';
  
  if (cds.options.Eigenvalues==1)
      varargout{2} = [varargout{2}; D];
  end

%-------------------------------------------------------
  
function option = options
global homds 
  % Check for symbolic derivatives in odefile
  
  symord = 0; 
  if ~isempty(homds.Jacobian), symord = 1; end
  if ~isempty(homds.Hessians), symord = 2; end
  if ~isempty(homds.Der3), symord = 3; end
  if  ~isempty(homds.Der4), symord = 4; end
  if ~isempty(homds.Der5), symord = 5; end

  option = contset;
  option = contset(option, 'SymDerivative', symord);
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', zeros(1,13));
  symordp = 0;
  if  ~isempty(homds.JacobianP), symordp = 1; end
  if ~isempty(homds.HessiansP),  symordp = 2; end
  option = contset(option, 'SymDerivativeP', symordp);
  
%------------------------------------------------------  
  
function [out, failed] = testf(id, x0, v)
global homds cds BigJac OldBigY

[x,x0,p,T,eps0,eps1,YS,YU] = rearr(x0);
ups = reshape(x,homds.nphase,homds.tps);

A = cjac(homds.func,homds.Jacobian,x0,num2cell(p),homds.ActiveParams);
D = eig(A);
nneg = sum(real(D) < 0);

[D1,indlist] = sort(real(D));
[val, ind] = min(abs(real(D1)));

if real(D1(ind)) > 0
    lambda1 = ind;
    lambda2 = ind+1;
    lambda3 = ind+2;
    mu1 = ind-1;
    mu2 = ind-2;
    mu3 = ind-3;
else
    lambda1 = ind+1;
    lambda2 = ind+2;
    lambda3 = ind+3;
    mu1 = ind;
    mu2 = ind-1;
    mu3 = ind-2;
end

D2 = D1;
D2(ind) = [];
[val2,ind2] = min(abs(real(D2)));
if ind2 >= ind
    ind2 = ind2 + 1;
end

res = zeros(14,1);
out = res;
failed = [];
% return

% at least 1 negative and 1 positive eigenvalue
if (mu1 > 0) && (lambda1 <= homds.nphase)
    % 1. Neutral saddle, saddle-focus or bi-focus
    res(1,1) = real(D1(mu1)) + real(D1(lambda1));

    % 2 negative eigenvales and 1 positive
    if (mu2 > 0)   
        % 2. Double real stable leading eigenvalues
        if imag(D1(mu1)) < cds.options.VarTolerance
            res(2,1) = (real(D1(mu1)) - real(D1(mu2)))^2;
        else
            res(2,1) = -(imag(D1(mu1)) - imag(D1(mu2)))^2;
        end
        
        % 4. Neutrally-divergent saddle-focus (stable)
        res(4,1) = real(mu1) + real(mu2) + real(lambda1);
        
        if (mu3 > 0)
            % 6. Three leading eigenvalues (stable)
            res(6,1) = real(mu1) - real(mu3);
        end
    end
    
    if (lambda2 <= homds.nphase)
        % 3. Double real unstable leading eigenvalues
        if imag(D1(lambda1)) < cds.options.VarTolerance
            res(3,1) = (real(D1(lambda1)) - real(D1(lambda2)))^2;
        else
            res(3,1) = -(imag(D1(lambda1)) - imag(D1(lambda2)))^2;
        end
        
        % 5. Neutrally-divergent saddle-focus (unstable)
        res(5,1) = real(mu1) + real(lambda2) + real(lambda1);
        
        if (lambda3 > 0)
            % 7. Three leading eigenvalues (unstable)
            res(7,1) = real(mu1) - real(mu3);
        end
    end
elseif mu1 == 0
    res(1,1) = D1(lambda1);
    res(4,1) = real(D1(lambda1));
else
    res(1,1) = D1(mu1);
    res(4,1) = real(D1(mu1));
end

% 8. Eigenvalue with zero real part => non-hyperbolic equilibrium
if nneg == homds.nneg
    % The signs of the eigenvalues are still all the same, if an
    % eigenvalue is small enough, we have a zero real part.    
    res(8,1) = sign(real(D1(ind)))*real(D1(ind)) - 10*cds.options.VarTolerance;
    
else
    % One of the negative eigenvalues has turned positive
    res(8,1) = real(D1(ind));
end

% 10. Non-central HSN: The eigenvalue with smallest real part has zero imaginary part
if abs(imag(D1(ind))) < cds.options.VarTolerance
    res(9,1) = res(8,1);
else
    res(9,1) = imag(D1(ind));
end

% 11. Bogdanov-Takens point: 2nd smallest eigenvalue must be zero
%prod = D1(mu1) * D1(lambda1);
%res(11,1) = sign(prod)*prod - cds.options.VarTolerance;
%res(11,1) = sign(real(D1(ind2)))*real(D1(ind2)) - cds.options.VarTolerance;
res(10,1) = norm(D1(ind2)) - cds.options.VarTolerance;

if (D1(1) > 0) || (D1(end) < 0)
    % This is a saddle-node, not a saddle!, so BT, NCH or something should
    % have been detected
    res(11:14,1) = 0;
    return;
end

V = [];
if (mu1 > 0)
if abs(imag(D1(mu1))) < cds.options.VarTolerance
    % 11. Orbit-flip with respect to stable manifold
    [V,D] = eig(A');
    % We still have the ordering of the eigenvalues of A (and thus A') from
    % before
    VN = V(:,indlist);
    w1s = VN(:,mu1);
    res(11,1) = exp(-D1(mu1) * homds.T) * (w1s' * (ups(:,end) - x0));
end
end

if (lambda1 <= homds.nphase)
if abs(imag(D1(lambda1))) < cds.options.VarTolerance
    if isempty(V)
        [V,D] = eig(A');
        VN = V(:,indlist);
    end
    w1u = VN(:,lambda1);
    % 12. Orbit-flip with respect to unstable manifold
    res(12,1) = exp(D1(lambda1) * homds.T) * (w1u' * (ups(:,1) - x0));
end

end

V = [];
if (mu1 > 0) & (~isempty(YS))
    if abs(imag(D1(mu1))) < cds.options.VarTolerance
        % 13. Inclination-flip with respect to stable manifold
        [V,D] = eig(A);
        VN = V(:,indlist);
        if homds.nneg == 0
            homds.nneg = homds.nneg + 1;
        elseif homds.nneg == homds.nphase
            homds.nneg = homds.nneg - 1;
        end
        
        BigJac2 = BigJac([1:homds.ncoords-homds.nphase end-1-homds.nphase:end-2],1:homds.ncoords);
        
        if (isempty(OldBigY)) | (size(OldBigY,2) < 2) | (size(OldBigY,1) ~= homds.ncoords)
            psiold = rand(homds.ncoords,1);
            phiold = rand(homds.ncoords,1);
            BigJacT = [BigJac2 psiold; phiold' 0];
            while condest(BigJacT) > 1e10
                psiold = rand(homds.ncoords,1);
                phiold = rand(homds.ncoords,1);
                BigJacT = [BigJac2 psiold; phiold' 0];
            end    
        else
            psiold = OldBigY(:,1);
            phiold = OldBigY(:,2);  
            BigJacT = [BigJac2 psiold; phiold' 0];         
        end
        
        phitemp = BigJacT \ [zeros(homds.ncoords,1);1];
        phinew = phitemp(1:end-1);
        
        psitemp = BigJacT' \ [zeros(homds.ncoords,1);1];
        psinew = psitemp(1:end-1);
        
        OldBigY = [psinew/norm(psinew) phinew/norm(phinew)];
        sigma2 = psinew(end-1);
        
        % Needed small eigenvector
        v1s = VN(:,mu1);
        % Orthogonal complement of stable eigenspace
        Q0S = homds.oldStableQ;
        LSorth = Q0S * [-YS'; eye(size(YS,1))];
        res(13,1) = exp(-D1(mu1) * T) * (v1s' * LSorth * sigma2);
        res(13,1) = sign(res(13,1)) * res(13,1) - cds.options.VarTolerance;
    end
else
    if lambda1 <= homds.nphase
        res(13,1) = -D1(lambda1);
        res(14,1) = D1(lambda1);
    else 
        res(13,1) = -D1(mu1);
        res(14,1) = D1(mu1);
    end
        
    failed = [];
    return
end

if (lambda1 <= homds.nphase) & (~isempty(YU))
    if abs(imag(D1(lambda1))) < cds.options.VarTolerance
        % 14. Inclination-flip with respect to unstable manifold
        if isempty(V)
            [V,D] = eig(A);
            VN = V(:,indlist);
            if homds.nneg == 0
                homds.nneg = homds.nneg + 1;
            elseif homds.nneg == homds.nphase
                homds.nneg = homds.nneg - 1;
            end

            if (isempty(OldBigY)) | (size(OldBigY,2) < 2) | (size(OldBigY,1) ~= homds.ncoords)
                psiold = rand(homds.ncoords,1);
                phiold = rand(homds.ncoords,1);
                BigJacT = [BigJac2 psiold; phiold' 0];
                while cond(BigJacT) > 1e10
                    psiold = rand(homds.ncoords,1);
                    phiold = rand(homds.ncoords,1);
                    BigJacT = [BigJac2 psiold; phiold' 0];
                end            
            else
                psiold = OldBigY(:,1);
                phiold = OldBigY(:,2);   
                BigJacT = [BigJac2 psiold; phiold' 0];         
            end
        
            phitemp = BigJacT \ [zeros(homds.ncoords,1);1];
            phinew = phitemp(1:end-1);
        
            psitemp = BigJacT' \ [zeros(homds.ncoords,1);1];
            psinew = psitemp(1:end-1);
        
            OldBigY = [psinew/norm(psinew) phinew/norm(phinew)];
        end
        
        sigma1 = psinew(end);
        
        % Needed small eigenvector
        v1u = VN(:,lambda1);
        % Orthogonal complement of stable eigenspace
        Q0U = homds.oldUnstableQ;
        LUorth = Q0U * [-YU'; eye(size(YU,1))];
        res(14,1) = exp(D1(lambda1) * T) * (v1u' * LUorth * sigma1);
        res(14,1) = sign(res(14,1)) * res(14,1) - cds.options.VarTolerance;
    end
else
    res(14,1) = D1(mu1);
end

out = res';
failed = [];

%-------------------------------------------------------------

function [out, failed] = userf(userinf, id, x, v)
global  homds
dim =size(id,2);
failed = [];
[x,x0,p,T,eps0,eps1,YS,YU] = rearr(x); p = num2cell(p);
out(dim) = 0;
for i=1:dim
  lastwarn('');
  if (userinf(i).state==1)
      out(i)=feval(homds.user{id(i)},0,x0,p{:});
  else
      out(i)=0;
  end
  if ~isempty(lastwarn)
    failed = [failed i];
  end
end

%-----------------------------------------------------------------

function [failed,s] = process(id, x, v, s)
global homds
[x,x0,p,T,eps0,eps1,YS,YU] = rearr(x);
AP = p(homds.ActiveParams);
failed = 0;
switch(id)    
    case 1
        A = cjac(homds.func,homds.Jacobian,x0,num2cell(p),homds.ActiveParams);
        D = eig(A);
        [D1,indlist] = sort(real(D));
        [val, ind] = min(abs(real(D1)));
    
        if real(D1(ind)) > 0
            lambda1 = ind;
            lambda2 = ind+1;
            lambda3 = ind+2;
            mu1 = ind-1;
            mu2 = ind-2;
            mu3 = ind-3;
        else
            lambda1 = ind+1;
            lambda2 = ind+2;
            lambda3 = ind+3;
            mu1 = ind;
            mu2 = ind-1;
            mu3 = ind-2;
        end
        if (imag(D1(mu1)) == 0) && (imag(D1(lambda1)) == 0)
            fprintf('Neutral saddle, parameters = %g and %g.\n',AP(1),AP(2));
            s.msg  = sprintf('Neutral saddle');
        elseif (imag(D1(mu1)) == 0) || (imag(D1(lambda1)) == 0)
            fprintf('Saddle-focus, parameters = %g and %g.\n',AP(1),AP(2));
            s.msg  = sprintf('Saddle-focus');
        else
            fprintf('Bi-focus, parameters = %g and %g.\n',AP(1),AP(2));
            s.msg  = sprintf('Bi-focus');
        end
    case 2
        fprintf('Double real stable leading eigenvalues, parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Double real stable leading eigenvalues (saddle - saddle-focus transition)');
    case 3
        fprintf('Double real unstable leading eigenvalues, parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Double real unstable leading eigenvalues (saddle - saddle-focus transition)');
    case 4
        fprintf('Neutrally-divergent saddle-focus (stable), parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Neutrally-divergent saddle-focus (stable)');
    case 5
        fprintf('Neutrally-divergent saddle-focus (unstable), parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Neutrally-divergent saddle-focus (unstable)');
    case 6
        fprintf('Three leading eigenvalues (stable), parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Three leading eigenvalues (stable)');
    case 7
        fprintf('Three leading eigenvalues (unstable), parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Three leading eigenvalues (unstable)');
    case 8
        fprintf('Shil''nikov-Hopf (Non-hyperbolic equilibrium), parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Shil''nikov-Hopf (Non-hyperbolic equilibrium)');
    case 9
        fprintf('Non-central homoclinic-to-saddle-node (Non-hyperbolic equilibrium), parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Non-central homoclinic-to-saddle-node (Non-hyperbolic equilibrium)');
    case 10
        fprintf('Bogdanov-Takens point, parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Bogdanov-Takens point');
    case 11
        fprintf('Orbit-flip with respect to the stable manifold, parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Orbit-flip with respect to the stable manifold');
    case 12
        fprintf('Orbit-flip with respect to the unstable manifold, parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Orbit-flip with respect to the unstable manifold');
    case 13
        fprintf('Inclination-flip with respect to the stable manifold, parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Inclination-flip with respect to the stable manifold');
    case 14
        fprintf('Inclination-flip with respect to the unstable manifold, parameters = %g and %g.\n',AP(1),AP(2));
        s.msg  = sprintf('Inclination-flip with respect to the unstable manifold');
end

%-------------------------------------------------------------  

function [S,L] = singmat

  S = 8 * (ones(14,14) - eye(14,14));
  S(1,8) = 1;
  S(8,[9 10]) = 1;
  S(9,8) = 0;
  S(9,10) = 1;
  L = ['NS ';'DRS';'DRU';'NDS';'NDU';'3LS';'3LU';'SH ';'NCH';'BT ';'OFS';'OFU';'IFS';'IFU'];

%--------------------------------------------------------

function locate(id, varargin)
error('No locator defined for singularity %d', id);
    
%----------------------------------------------------------

function varargout = init(varargin)

WorkspaceInit(varargin{1:2});
% all done succesfully
varargout{1} = 0;

%-----------------------------------------------------------

function done

%-----------------------------------------------------------

function [res,x,v] = adapt(x,v)
global homds cds

cds.adapted = 1;

YU = homds.YU;
YS = homds.YS;

[x,v] = Hom_adapt_mesh(x,v);

Q0S = homds.oldStableQ;
QbS1 = Q0S(:,1:homds.nneg);
QbS2 = Q0S(:,homds.nneg+1:end);
Q0U = homds.oldUnstableQ;
QbU1 = Q0U(:,1:homds.npos);
QbU2 = Q0U(:,homds.npos+1:end);

[Q1S,S1,R1] = svd(QbS1 + QbS2*YS);
[Q1U,S1,R1] = svd(QbU1 + QbU2*YU);

Q1S = Q1S*diag(sign(diag(homds.oldStableQ'*Q1S)));
Q1U = Q1U*diag(sign(diag(homds.oldUnstableQ'*Q1U)));

homds.oldStableQ = Q1S;
homds.oldUnstableQ = Q1U;

homds.YS = zeros(size(homds.YS));
homds.YU = zeros(size(homds.YU));

x(end-2*homds.nneg*homds.npos+1:end)=zeros(2*homds.nneg*homds.npos,1);
%v(end-2*homds.nneg*homds.npos+1:end)=zeros(2*homds.nneg*homds.npos,1);
res = 1;

%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------
 
function [x,x0,p,T,eps0,eps1,YS,YU] = rearr(x1)
% Rearranges x1 into all of its components
global homds

x = x1(homds.coords);
x0 = x1(homds.ncoords+1:homds.ncoords+homds.nphase);

p = homds.P0;
p(homds.ActiveParams) = x1(homds.PeriodIdx+1:homds.PeriodIdx+2);
idx = homds.PeriodIdx+3;

if homds.extravec(1)
    T = x1(idx);
    idx = idx+1;
else
    T = homds.T;
end
if homds.extravec(2)
    eps0 = x1(idx);
    idx = idx+1;
else
    eps0 = homds.eps0;
end
if homds.extravec(3)
    eps1 = x1(idx);
    idx = idx+1;
else
    eps1 = homds.eps1;
end

    YU = reshape(x1(idx:idx+homds.npos*homds.nneg-1),homds.nneg,homds.npos);
    idx = idx + homds.npos*homds.nneg;
    YS = reshape(x1(idx:idx+homds.npos*homds.nneg-1),homds.npos,homds.nneg);
    
% -------------------------------------------------------------

% ---------------------------------------------------------------

function WorkspaceInit(x,v)
global homds
homds.cols_p1 = 1:(homds.ncol+1);
homds.cols_p1_coords = 1:(homds.ncol+1)*homds.nphase;
homds.ncol_coord = homds.ncol*homds.nphase;
homds.col_coords = 1:homds.ncol*homds.nphase;
homds.coords = 1:homds.ncoords;
homds.pars = homds.ncoords+(1:2);
homds.tsts = 1:homds.ntst;
homds.cols = 1:homds.ncol;
homds.phases = 1:homds.nphase;
homds.ntstcol = homds.ntst*homds.ncol;

homds.idxmat = reshape(fix((1:((homds.ncol+1)*homds.ntst))/(1+1/homds.ncol))+1,homds.ncol+1,homds.ntst);
homds.dt = homds.msh(homds.tsts+1)-homds.msh(homds.tsts);

homds.wp = kron(homds.wpvec',eye(homds.nphase));
homds.pwwt = kron(homds.wt',eye(homds.nphase));
homds.pwi = homds.wi(ones(1,homds.nphase),:);

homds.wi = nc_weight(homds.ncol)';

[homds.bialt_M1,homds.bialt_M2,homds.bialt_M3,homds.bialt_M4]=bialtaa(homds.nphase);

% ------------------------------------------------------

function [x,v,s] = WorkspaceDone(x,v,s)

%------------------------------------------------------------
