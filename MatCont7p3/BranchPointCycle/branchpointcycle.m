function out = branchpointcycle
%
% BPCycles curve definition file
% 

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

%--------------------------------------------------------
function func = curve_func(arg)
global lds 
[x,p,T] = rearr(arg);
  f = BVP('BVP_LC_f','BVP_LC_bc','BVP_LC_ic',x,p,T);
  % append g
  J = BVP_BPCjac('BVP_BPC_jacCC',x,p,T,1,1);
  
  b=[zeros(lds.ncoords+1,2);eye(2)];
  s = J\b;
  if lds.BPC_switch==1
      lds.BPC_new_phi = s(1:lds.ncoords+2,:)';
  else
      b = []; b(lds.ncoords+3)=1;
      t=J'\b';
      lds.BPC_new_psi=t(1:lds.ncoords+1)';
  end
  f(end+1:end+2) = s(end,:);
  func = f;
%----------------------------------------------------------  
function varargout = jacobian(varargin)
global  lds
  [x,p,T] = rearr(varargin{1});
  f = BVP_jac('BVP_LCX_jac',x,p,T,4,2);
  
  % calculate v and w
  j = size(f,1)+1;
  f(j:j+1,:) = 0;
  b=[zeros(lds.ncoords+1,2);eye(2)];
  J = BVP_BPCjac('BVP_BPC_jacCC',x,p,T,1,1);
  
  sn = J\b;
  b = []; b(lds.ncoords+3)=1;
  st = J'\b';
  v11 = sn(1:lds.ncoords,1)';
  v21 = sn(1:lds.ncoords,2)';
  v12 = sn(lds.ncoords+1,1);
  v22 = sn(lds.ncoords+1,2);
  v13 = sn(lds.ncoords+2,1);
  v23 = sn(lds.ncoords+2,2);
  w1 = st(1:lds.ncoords-lds.nphase)';
  % calculate g'
  ups = reshape(x,lds.nphase,lds.tps);
  p = num2cell(p);

  
  range0 = lds.cols_p1;
  range1 = lds.col_coords;
  range2 = lds.cols_p1_coords;
  
  t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
  kr1 = fix(t/lds.nphase);
  kr2 = rem(t,lds.nphase)+1;
  for tstpt=lds.tsts
    xp  = ups(:,range0)*lds.wt;
    cv1 = v11(range2)';
    cv2 = v21(range2)';
    cw1 = w1(range1);
   
    range = lds.phases;
    for c=lds.cols
      xt = xp(:,c);
      sysj   = cjac(lds.func,lds.Jacobian,xt,p,lds.ActiveParams);
      sysp   = cjacp(lds.func,lds.JacobianP,xt,p,lds.ActiveParams);
      sysh   = chess(lds.func,lds.Jacobian,lds.Hessians,xt,p,lds.ActiveParams);
      syshp  = chessp(lds.func,lds.Jacobian,lds.HessiansP,xt,p,lds.ActiveParams);
      sysbr  = cjacbr(lds.func,lds.JacobianP,xt,p,lds.ActiveParams,lds.BranchParam);
      syshbr = chessbr(lds.func,lds.Jacobian,lds.HessiansP,xt,p,lds.ActiveParams,lds.BranchParam);
      syshbrp = chesspbr(lds.func,lds.JacobianP,lds.HessiansP,xt,p,lds.ActiveParams,lds.BranchParam);
      wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
      for d=lds.phases
        sh1(:,d) = (wtk.*sysh(:,kr2,d))*cv1;
        sh2(:,d) = (wtk.*sysh(:,kr2,d))*cv2;
      end      
      t11 = T* wtk.*sh1(:,kr2) + wtk.*sysj(:,kr2)*v12 + T*wtk.*syshbr(:,kr2,1)*v13;
      t12 = T* wtk.*sh2(:,kr2) + wtk.*sysj(:,kr2)*v22 + T*wtk.*syshbr(:,kr2,1)*v23;
      t21 = (wtk.*sysj(:,kr2))*cv1 + sysbr*v13;
      t22 = (wtk.*sysj(:,kr2))*cv2 + sysbr*v23;
      t31 = T* wtk.*syshp(:,kr2,1)* cv1 + sysp(:,1)*v12 + T*syshbrp(:,:,1)*v13;
      t32 = T* wtk.*syshp(:,kr2,1)* cv2 + sysp(:,1)*v22 + T*syshbrp(:,:,1)*v23;
      t41 = T* wtk.*syshp(:,kr2,2)* cv1 + sysp(:,2)*v12 + T*syshbrp(:,:,2)*v13;
      t42 = T* wtk.*syshp(:,kr2,2)* cv2 + sysp(:,2)*v22 + T*syshbrp(:,:,2)*v23;
      t51 = T* wtk.*syshp(:,kr2,3)* cv1 + sysp(:,3)*v12 + T*syshbrp(:,:,3)*v13;
      t52 = T* wtk.*syshp(:,kr2,3)* cv2 + sysp(:,3)*v22 + T*syshbrp(:,:,3)*v23;

      syshess1(range,:) = [t11 t21 t31 t41 t51];      
      syshess2(range,:) = [t12 t22 t32 t42 t52];
      range = range + lds.nphase;
    end
    f(j,[range2 lds.ncoords+(1:4)])   = f(j,[range2 lds.ncoords+(1:4)])   + cw1*syshess1;
    f(j+1,[range2 lds.ncoords+(1:4)]) = f(j+1,[range2 lds.ncoords+(1:4)]) + cw1*syshess2;
    range0 = range0 + lds.ncol;
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol_coord;
  end
  varargout{1} = f;
%-------------------------------------------------------
function hessians(varargin)  
%-------------------------------------------------------
function varargout = defaultprocessor(varargin)
global lds cds
  [x,p,T] = rearr(varargin{1});
  v = rearr(varargin{2});

  % update
  lds.ups = reshape(x,lds.nphase,lds.tps);
  lds.vps = reshape(v,lds.nphase,lds.tps);

  if nargin > 2
    % set data in special point structure
    s = varargin{3};
    s.data.timemesh = lds.msh;
    s.data.ntst = lds.ntst;
    s.data.ncol = lds.ncol;
    s.data.parametervalues = p;
    s.data.T = T;
    varargout{3} = s;
  end

  % update upoldp
  p = num2cell(p);
  for i=1:lds.tps
    lds.upoldp(:,i) = T*feval(lds.func, 0, lds.ups(:,i), p{:});
  end
  
  % calculate multipliers if requested
  if lds.CalcMultipliers %& (isempty(lds.multipliersX)|lds.multipliersX~=varargin{1})
    lds.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,varargin{1},[]));
    lds.multipliersX = varargin{1};
  end
  if lds.CalcMultipliers==0
      lds.multipliers = [];
  end
  % special data
  varargout{2} = [lds.msh'; lds.multipliers];

  % all done succesfully
  varargout{1} = 0;
%-------------------------------------------------------
function option = options(varargin)
global lds
  % Check for symbolic derivatives in odefile  
  option = contset;
  
  symord = 0; 
  if ~isempty(lds.Jacobian), symord = 1; end
  if ~isempty(lds.Hessians), symord = 2; end
  if ~isempty(lds.Der3), symord = 3; end
  if ~isempty(lds.Der4), symord = 4; end
  if ~isempty(lds.Der5), symord = 5; end
  option = contset(option, 'SymDerivative', symord);
  
  symordp = 0;
  if ~isempty(lds.JacobianP), symordp = 1; end
  if ~isempty(lds.HessiansP),  symordp = 2; end
  option = contset(option, 'SymDerivativeP', symordp);
  
  option = contset(option, 'Singularities', 0);
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', [0 0 0]);
  
%------------------------------------------------------------------
function testf(varargin)
%-------------------------------------------------------------------
function [out, failed] = userf(userinf, id, x, v)
global lds
dim =size(id,2);
failed = [];
out(dim) = 0;
for i=1:dim
  lastwarn('');
  [x0,p] = rearr(x); p = num2cell(p);
  if (userinf(i).state==1)
      out(i)=feval(lds.user{id(i)},0,x0,p{:});
  else
      out(i)=0;
  end
  if ~isempty(lastwarn)
    failed = [failed i];
  end
end
%------------------------------------------------------------
function process(varargin)  
%------------------------------------------------------------
function singmat(varargin)
%------------------------------------------------------------
function locate(varargin)
%------------------------------------------------------------
function varargout = init(varargin)
  WorkspaceInit(varargin{1:2});
  % all done succesfully
  varargout{1} = 0;
%-------------------------------------------------------------
function done
  WorkspaceDone;
%-------------------------------------------------------------
function [res,x,v] = adapt(x,v)
global lds
% calculate phi and psi for next point
if lds.BPC_switch == 0
  lds.BPC_phi = lds.BPC_new_phi/norm(lds.BPC_new_phi);
  lds.BPC_phi1 = lds.BPC_new_phi(1,:);
  lds.BPC_phi2 = lds.BPC_new_phi(2,:);
else
  lds.BPC_psi = lds.BPC_new_psi/norm(lds.BPC_new_psi);
end
lds.BPC_switch = 1-lds.BPC_switch;
[x,v]=adapt_mesh(x,v);
res = [];


%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------

function [x,p,T] = rearr(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x), parameters (p) and period (T)
global lds
p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+(1:3));
x = x0(lds.coords);
T = x0(lds.PeriodIdx);

% -------------------------------------------------------------

function f = BVP(BVP_f,BVP_bc,BVP_ic,x,p,T)
global lds

% extract ups
ups = reshape(x,lds.nphase,lds.tps);
p = num2cell(p);
f = zeros(1,lds.ncoords+1);
% function
range1 = lds.cols_p1;
range2 = lds.phases;
for j=lds.tsts
  xp = ups(:,range1)*lds.wt;
  t  = ups(:,range1)*lds.wpvec/lds.dt(j);
  for c=lds.cols
    f(range2) = feval(BVP_f,lds.func,t(:,c),xp(:,c),p,T);
    range2 = range2+lds.nphase;
  end
  range1 = range1+lds.ncol;
end

% boundary conditions
f(range2) = feval(BVP_bc,ups(:,1),ups(:,lds.tps));
% integral constraint
f(lds.ncoords+1) = feval(BVP_ic,ups);

f = f';
%------------------------------------------------------------------
function jac = BVP_BPCjac(BVP_func,x,p,T,pars,nc)
global lds
p2 = num2cell(p);
jac = feval(BVP_func,lds.func,x,p,T,pars,nc,lds,p2,lds.Jacobian,lds.ActiveParams,lds.JacobianP,lds.BranchParam); 

% -------------------------------------------------------------

function jac = BVP_jac(BVP_func,x,p,T,pars,nc)
global lds 
 
p2 = num2cell(p);
jac = feval(BVP_func,lds.func,x,p,T,pars,nc,lds,p2,lds.Jacobian,lds.ActiveParams,lds.JacobianP); 
jac = sparse(full(jac));


% ---------------------------------------------------------------


function WorkspaceInit(x,v)
global cds lds


lds.BPC_new_phi = [lds.BPC_phi1;lds.BPC_phi2];
lds.BPC_new_psi = lds.BPC_psi;
lds.BPC_switch = 0;

lds.CalcMultipliers = contget(cds.options, 'Multipliers', 0);
lds.multipliersX = [];
lds.multipliers = nan;

r = (0:(lds.ntst*lds.nphase-1));
lds.multi_r1 = (floor(r./lds.nphase)+1)*lds.ncol_coord-lds.nphase+mod(r,lds.nphase)+1;
r = (0:((lds.ntst+1)*lds.nphase-1));
lds.multi_r2 = floor(r./lds.nphase)*lds.ncol_coord+mod(r,lds.nphase)+1;
lds.pars = lds.ncoords+(1:4);
% ------------------------------------------------------

function WorkspaceDone
