function out = fixedperiodlimitcycle
%
% Limit cycle curve definition file for a problem in odefile
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

% lds.nphase = dimension of a single point on the cycle
% lds.ncol = number of collocation points
% lds.ntst = number of test points
% lds.tps = ncol*ntst+1


%----------------------------------------------------
function func = curve_func(arg)

  [x,p,T] = rearr(arg);
  func = BVP('BVP_LC_f','BVP_LC_bc','BVP_LC_ic',x,p,T);

%------------------------------------------------------
function varargout = jacobian(varargin)

  [x,p,T] = rearr(varargin{1});
  varargout{1} = BVP_jac('BVP_LCT_jac_f','BVP_LC_jac_bc','BVP_LC_jac_ic',x,p,T,2,2);
  
%-----------------------------------------------------

function hessians(varargin)
%------------------------------------------------------

function varargout = defaultprocessor(varargin)
global lds cds
  [x,p,T] = rearr(varargin{1});
  v = rearr(varargin{2});
  % update
  lds.ups = reshape(x,lds.nphase,lds.tps);
  lds.vps = reshape(v,lds.nphase,lds.tps);
  % update upoldp
  p1 = num2cell(p);
  for i=1:lds.tps
    lds.upoldp(:,i) = T*feval(lds.func, 0, lds.ups(:,i), p1{:});
  end
  % calculate multipliers if requested
  if lds.CalcMultipliers && (isempty(lds.multipliersX)||lds.multipliersX~=varargin{1})
      try
          lds.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,varargin{1},[])); 
          lds.multipliersX = varargin{1};
      catch
          lds.multipliersX =varargin{1};
      end
  end
  
  if nargin > 2
    % set data in special point structure
    s = varargin{3};
    s.data.multipliers = lds.multipliers;
    s.data.timemesh = lds.msh;
    s.data.ntst = lds.ntst;
    s.data.ncol = lds.ncol;
    s.data.parametervalues = p;
    s.data.T = T;
    s.data.phi = lds.PD_phi(lds.coords);
    varargout{3} = s;
  end
  if lds.CalcMultipliers==0
      lds.multipliers=nan;
  end
  % special data
  varargout{2} = [lds.multipliers;lds.PD_new_phi(:)];
  % all done succesfully
  varargout{1} = 0;
%-------------------------------------------------------
function option = options
global lds 
  % Check for symbolic derivatives in odefile
  
  symord = 0; 
  if ~isempty(lds.Jacobian), symord = 1; end
  if ~isempty(lds.Hessians), symord = 2; end
  if ~isempty(lds.Der3), symord = 3; end
  if ~isempty(lds.Der4), symord = 4; end
  if ~isempty(lds.Der5), symord = 5; end
  
  option = contset;
  switch lds.nphase
      case 1
          option=contset(option,'IgnoreSingularity',[2 3 4]);
      case 2
          option=contset(option,'IgnoreSingularity',[4]);
  end

  option = contset(option, 'SymDerivative', symord);
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', [1 0 0 0]);
  symordp = 0;
  if ~isempty(lds.JacobianP), symordp = 1; end
  if ~isempty(lds.HessiansP),  symordp = 2; end
  option = contset(option, 'SymDerivativeP', symordp);
  
%------------------------------------------------------  
function [out, failed] = testf(id, x0, v)
global lds cds

[x,p,T] = rearr(x0);
out(3)=0;
failed = [];
if any(ismember([6 8],id)) && (isempty(lds.multipliersX)||(lds.multipliersX~=x0))
    lds.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,x0,[])); 
    lds.multipliersX = x0;
end

  
for i=id
    lastwarn('');
    switch i
        
    case {1,2,3,4,5} %BPC
        Jb = BVP_jac('BVP_BPC_jac_f','BVP_BPC_jac_bc','BVP_BPC_jac_ic',x,p,T,1,1);
        Jb= [Jb [lds.BP_psi1(lds.coords)'; 0];lds.BP_phi1(lds.coords) 0 0];
        bb = [zeros(lds.ncoords,2);eye(2)] ;
        if lds.BP_switch == 1
             sp = Jb'\bb;
             lds.BP_new_psi = reshape(sp(lds.coords,1),lds.nphase,lds.tps);
             lds.BP_new_psi1= reshape(sp(lds.coords,2),lds.nphase,lds.tps);
        end
        
        st = Jb\bb;
        lds.BP_new_phi = reshape(st(lds.coords,1),lds.nphase,lds.tps);
        lds.BP_new_phi1= reshape(st(lds.coords,2),lds.nphase,lds.tps);        
        out(1) = st(end-1,1);
        out(2) = st(end-1,2);
        out(3) = st(end,1);
        out(4) = st(end,2);
        out(5) = norm(out(1:4));
    case 6 % PD
        A = lds.monodromy;
        A = A + eye(size(A,1));
        out(6) = det(A);
    case 7 % LPC
        out(7) = v(lds.ncoords+2);
    case 8 %NS
        if ~lds.CalcMultipliers |isempty(lds.multipliersX)|(lds.multipliersX~=x0)
            lds.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,x0,[])); 
            lds.multipliersX = x0;
        end
        A = lds.monodromy;
        A = A(lds.bialt_M1).*A(lds.bialt_M2)-A(lds.bialt_M3).*A(lds.bialt_M4);
        A = A-eye(size(A,1));   
        out(8) = det(A);

    otherwise
        error('No such testfunction');
    end
  if ~isempty(lastwarn)
    failed = [failed i];
  end

end
%-------------------------------------------------------------
function [out, failed] = userf( userinf, id, x, v)
global  lds
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
%-----------------------------------------------------------------
function [failed,s] = process(id, x, v, s)
global lds 
  switch id
  case 1
    fprintf('Branch Point cycle(parameters = [%e;%e])\n',x(length(x)-1),x(length(x)));
    s.msg  = sprintf('Branch Point cycle'); 
  case 2
    [x0,p,T] = rearr(x);  
    J = BVP_jac('BVP_PD_jac_f','BVP_PD_jac_bc','BVP_PD_jac_ic',x0,p,T,1,1);
    [LJ,UJ] = lu(J);
    b = []; b(lds.ncoords+1)=1; b=b';
    ss = UJ\(LJ\b);
    lds.PD_phi = reshape(ss(lds.coords),lds.nphase,lds.tps);
    s.data.phi = lds.PD_phi(lds.coords);
    s.data.pdcoefficient = nf_PD(x);
    fprintf('Period Doubling (parameters = [%e;%e])\n',x(length(x)-1),x(length(x)));
    s.msg  = sprintf('Period Doubling');
    fprintf('Normal form coefficient = %d\n', s.data.pdcoefficient);
  case 3
    s.data.lpccoefficient = nf_LPC(x);
    fprintf('Limit point cycle (parameters = [%e;%e])\n',x(length(x)-1),x(length(x)));
    s.msg  = sprintf('Limit point cycle'); 
    fprintf('Normal form coefficient = %d\n', s.data.lpccoefficient);
  case 4
    s.data.nscoefficient = nscoefficient(x);
    fprintf('Neimark-Sacker (parameters = [%e;%e])\n',x(length(x)-1),x(length(x)));
    s.msg  = sprintf('Neimark Sacker');  
    fprintf('Normal form coefficient = %d\n', s.data.nscoefficient);
  end
  failed = 0;
%-------------------------------------------------------------  
function [S,L] = singmat

  S = [ 0 0 0 0 8 8 8 8
        8 8 8 8 8 0 8 8
        8 8 8 8 8 8 0 8
        8 8 8 8 1 8 1 0];


  L = [ 'BPC';'PD '; 'LPC'; 'NS ' ];

 
%--------------------------------------------------------
function [x,v] = locate(id, x1, v1, x2, v2)
switch id   
  case 1
    [x,v] = locateBPC(id, x1, v1, x2, v2);
  otherwise
    error('No locator defined for singularity %d', id);
end
%----------------------------------------------------------
function varargout = init(varargin)

  WorkspaceInit(varargin{1:2});
  % all done succesfully
  varargout{1} = 0;
%-----------------------------------------------------------
function done

%-----------------------------------------------------------
function [res,x,v] = adapt(x,v)
global lds 

% calculate phi and psi for next point
if lds.BP_switch == 0
  lds.BP_phi = lds.BP_new_phi;
  lds.BP_phi =  lds.BP_phi/norm(lds.BP_phi(lds.coords));
  lds.BP_phi1 = lds.BP_new_phi1;
  lds.BP_phi1 = (lds.BP_phi1*lds.BP_phi')*lds.BP_phi-(lds.BP_phi*lds.BP_phi')*lds.BP_phi1;
  lds.BP_phi1 = lds.BP_phi1/norm(lds.BP_phi1(lds.coords));  
else
  lds.BP_psi = lds.BP_new_psi;
  lds.BP_psi =  lds.BP_psi/norm(lds.BP_psi);
  lds.BP_psi1 = lds.BP_new_psi1;
  lds.BP_psi1 = (lds.BP_psi1*lds.BP_psi')*lds.BP_psi-(lds.BP_psi*lds.BP_psi')*lds.BP_psi1;
  lds.BP_psi1 = lds.BP_psi1/norm(lds.BP_psi1);
end
lds.BP_switch = 1-lds.BP_switch;

[x,v] = adapt_mesh(x,v);
res = 1;



%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------


function [x,v] = locateBPC(id, x1, v1, x2, v2)
global  cds

ndim = cds.ndim;

initpq(x1);
b = 0;
x = x1;
i = 0;

v = 0.5*(v1+v2);
u = [x; b];
[A,f]=locjac(x,b);
while i < 4
  du = A\f;
  u = u - du;

  x = u(1:ndim);
  b = u(ndim+1);

  [A,f]=locjac(x,b);
  % WM: VarTol and FunTol were switched
  if norm(du) < cds.options.VarTolerance && norm(f) < cds.options.FunTolerance
      return; 
  end
  
  i = i+1;
     
end
x = [];


% ---------------------------------------------------------------

function [A, f] = locjac(x0, b)
% A = jac of system
% f = system evaluated at (x,b)
global cds lds

[x,p,T] = rearr(x0);

% append g
J = BVP_BPC_jac('BVP_BPCT_f','BVP_BPC_bc1',x,p,T,1,1);
b1=[zeros(lds.ncoords+1,2);eye(2)];
sn = J\b1;
f = [feval(cds.curve_func, x0) + b*lds.BPC_psi'; sn(end,:)'];

A = BVP_jac('BVP_LCT_jac_f','BVP_LC_jac_bc','BVP_LC_jac_ic',x,p,T,2,2);

j = size(A,1)+1;
A(:,j+1) = lds.BPC_psi';
A(j,:)   = 0;
A(j+1,:) = 0;

b1 = []; b1(lds.ncoords+3)=1;
st = J'\b1';
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

for tstpt = lds.tsts
    xp  = ups(:,range0)*lds.wt;
    cv1 = v11(range2)';
    cv2 = v21(range2)';
    cw1 = w1(range1);
    
    range = lds.phases;
    for c=lds.cols
        xt = xp(:,c);

        sysh   = chess(lds.func,lds.Jacobian,lds.Hessians,xt,p,lds.ActiveParams);
        syshp  = chessp(lds.func,lds.Jacobian,lds.HessiansP,xt,p,lds.ActiveParams);
        syshbr = chessbr(lds.func,lds.Jacobian,lds.HessiansP,xt,p,lds.ActiveParams,lds.BranchParam);
        syshbrp = chesspbr(lds.func,lds.JacobianP,lds.HessiansP,xt,p,lds.ActiveParams,lds.BranchParam);
        wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
        for d=lds.phases
            sh1(:,d) = (wtk.*sysh(:,kr2,d))*cv1;
            sh2(:,d) = (wtk.*sysh(:,kr2,d))*cv2;
        end      
        t11 = T * wtk.*sh1(:,kr2) + T*wtk.*syshbr(:,kr2,1)*v12 + T*wtk.*syshbr(:,kr2,2)*v13;
        t12 = T* wtk.*sh2(:,kr2) + T*wtk.*syshbr(:,kr2,1)*v22 + T*wtk.*syshbr(:,kr2,2)*v23;
        t21 = T* wtk.*syshp(:,kr2,1)* cv1 + T*syshbrp(:,1,1)*v12 + T*syshbrp(:,1,2)*v13;
        t22 = T* wtk.*syshp(:,kr2,1)* cv2 + T*syshbrp(:,1,1)*v22 + T*syshbrp(:,1,2)*v23;
        t31 = T* wtk.*syshp(:,kr2,2)* cv1 + T*syshbrp(:,2,1)*v12 + T*syshbrp(:,2,2)*v13;
        t32 = T* wtk.*syshp(:,kr2,2)* cv2 + T*syshbrp(:,2,1)*v22 + T*syshbrp(:,2,2)*v23;
        syshess1(range,:) = [t11 t21 t31 zeros(size(t11,1),1)];      
        syshess2(range,:) = [t12 t22 t32 zeros(size(t12,1),1)];
        range = range + lds.nphase;
    end
    A(j,[range2 lds.ncoords+(1:3)])   = A(j,[range2 lds.ncoords+(1:3)])   + cw1*syshess1;
    A(j+1,[range2 lds.ncoords+(1:3)]) = A(j+1,[range2 lds.ncoords+(1:3)]) + cw1*syshess2;
    range0 = range0 + lds.ncol;
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol_coord;
end

%--------------------------------------------------------------
function  initpq(x0)
% A = jac of system
% f = system evaluated at (x,b,p)
global lds

[x,p,T] = rearr(x0);

% append g
ups = reshape(x,lds.nphase,lds.tps);
p = num2cell(p);
pars1 = lds.ncoords+1;
pars2 = lds.ncoords+2;

jac = spalloc(lds.ncoords+1,lds.ncoords+2,(lds.ncol+4)*lds.nphase);
% function
range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
for j = lds.tsts
  xp = ups(:,range0)*lds.wt;
  jac(range1,[range2 pars1 pars2]) = bordBVP_BPCT_f(lds.odefile,xp,p,T,j);
  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end
% boundary conditions
range  = (lds.tps-1)*lds.nphase+ (lds.phases);
range1 = lds.ncoords-lds.nphase+lds.phases;
jac(range,[lds.phases range1]) = bordBVP_LPC_bc1;
% integral constraint
ic = zeros(1,lds.ncoords);
range1 = lds.cols_p1;
range2 = lds.cols_p1_coords;
for j=lds.tsts
  pt = lds.dt(j)*(ups(:,range1).*lds.pwi);
  ic(range2) = ic(range2)+pt(lds.cols_p1_coords);
  range1 = range1 + lds.ncol;
  range2 = range2 + lds.ncol_coord;
end
jac(range(end)+1,1:lds.ncoords)= ic;
%compute borders
[Q,R,E] = qr(full(jac));
R(end,end) = 0;R(end,end-1) = 0;
p = E*[R(1:end-1,1:end-2)\-R(1:end-1,end-1:end);eye(2)];
lds.BPC_phi = p'/norm(p);
p = Q(:,end);
lds.BPC_psi = p';



%--------------------------------------------------------------
function [x,p,T] = rearr(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x) and parameters (p)
global lds

p = lds.P0;
T = lds.T;
p(lds.ActiveParams) = x0(end-1:end);
x = x0(lds.coords);

% -------------------------------------------------------------

function f = BVP(BVP_f,BVP_bc,BVP_ic,x,p,T)
global lds 

% extract ups
ups = reshape(x,lds.nphase,lds.tps);
p = num2cell(p);

% function
range1 = lds.cols_p1;
range2 = lds.phases;
for j=lds.tsts
  % value of polynomial on each collocation point
  xp = ups(:,range1)*lds.wt;
  
  % derivative of polynomial on each collocation point
  t  = ups(:,range1)*lds.wpvec/lds.dt(j);

  % evaluate function value on each collocation point
  for c=lds.cols
    f(range2) = feval(BVP_f,lds.func,t(:,c),xp(:,c),p,T);
    range2 = range2+lds.nphase;
  end

  range1 = range1+lds.ncol;
end
% boundary conditions
f(range2) = feval(BVP_bc,ups(:,1),ups(:,lds.tps));

% integral constraint
f(end+1) = feval(BVP_ic,ups);

f = f';

% -------------------------------------------------------------

function jac = BVP_jac(BVP_jac_f,BVP_jac_bc,BVP_jac_ic,x,p,T,pars,nc)
global lds 

ups = reshape(x,lds.nphase,lds.tps);
p = num2cell(p);
pars = lds.ncoords+(1:pars);

jac = spalloc(lds.ncoords+1,lds.ncoords+1,(lds.ncol+4)*lds.nphase);
%jac = zeros(lds.ncoords+1,lds.ncoords+length(p)-1);

range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;

for j=lds.tsts
  % value of polynomial on each collocation point
  xp = ups(:,range0)*lds.wt;

  % evaluate part of Jacobian matrix
  jac(range1,[range2 pars]) = feval(BVP_jac_f,lds.func,xp,p,T,j);

  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end

% boundary conditions
range = (lds.tps-1)*lds.nphase + (lds.phases);
jac(range,[lds.phases range lds.ncoords+(1:nc)]) = feval(BVP_jac_bc);

% integral constraint
jac(end,[lds.coords]) = feval(BVP_jac_ic);

% ---------------------------------------------------------------
function WorkspaceInit(x,v)
global cds lds
lds.cols_p1 = 1:(lds.ncol+1);
lds.cols_p1_coords = 1:(lds.ncol+1)*lds.nphase;
lds.ncol_coord = lds.ncol*lds.nphase;
lds.col_coords = 1:lds.ncol*lds.nphase;
lds.coords = 1:lds.ncoords;
lds.pars = lds.ncoords+(1:2);
lds.tsts = 1:lds.ntst;
lds.cols = 1:lds.ncol;
lds.phases = 1:lds.nphase;
lds.ntstcol = lds.ntst*lds.ncol;

lds.idxmat = reshape(fix((1:((lds.ncol+1)*lds.ntst))/(1+1/lds.ncol))+1,lds.ncol+1,lds.ntst);
lds.dt = lds.msh(lds.tsts+1)-lds.msh(lds.tsts);

lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwwt = kron(lds.wt',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);

lds.wi = nc_weight(lds.ncol)';


lds.PD_psi = reshape(exp((1:lds.ncoords)/lds.ncoords),lds.nphase,lds.tps);
lds.PD_psi = lds.PD_psi/norm(lds.PD_psi(lds.coords));
lds.PD_phi = reshape(ones(lds.ncoords,1),lds.nphase,lds.tps);
lds.PD_phi = lds.PD_phi/sqrt(BVP_PD_jac_ic*lds.PD_phi(lds.coords)');

lds.PD_new_phi = lds.PD_phi;
lds.PD_new_psi = lds.PD_psi;
lds.PD_switch = 0;

lds.BP_psi = reshape(exp((1:lds.ncoords)/lds.ncoords),lds.nphase,lds.tps);
lds.BP_psi = lds.BP_psi/norm(lds.BP_psi(lds.coords));
lds.BP_phi = reshape(ones(lds.ncoords,1),lds.nphase,lds.tps);
lds.BP_phi = lds.BP_phi/sqrt(BVP_BPC_jac_ic*lds.BP_phi(lds.coords)');

lds.BP_psi1 = reshape(ones(lds.ncoords,1),lds.nphase,lds.tps);
lds.BP_psi1 = (lds.BP_psi1*lds.BP_psi')*lds.BP_psi-(lds.BP_psi*lds.BP_psi')*lds.BP_psi1;
lds.BP_psi1 = lds.BP_psi1/norm(lds.BP_psi1(lds.coords));

lds.BP_phi1 = reshape(exp((1:lds.ncoords)/lds.ncoords),lds.nphase,lds.tps);
lds.BP_phi1 = (lds.BP_phi1*lds.BP_phi')*lds.BP_phi-(lds.BP_phi*lds.BP_phi')*lds.BP_phi1;
lds.BP_phi1 = lds.BP_phi1/norm(lds.BP_phi1(lds.coords));

lds.BP_new_phi = lds.BP_phi;
lds.BP_new_psi = lds.BP_psi;
lds.BP_new_psi1 = lds.BP_psi1;
lds.BP_new_phi1 = lds.BP_phi1;
lds.BP_switch = 0;

lds.BPC_switch = 0;

[lds.bialt_M1,lds.bialt_M2,lds.bialt_M3,lds.bialt_M4]=bialtaa(lds.nphase);
lds.CalcMultipliers = contget(cds.options, 'Multipliers', 0);
lds.multipliersX = [];
lds.multipliers = nan;
lds.monodromy = [];

r = (0:(lds.ntst*lds.nphase-1));
lds.multi_r1 = (floor(r./lds.nphase)+1)*lds.ncol_coord-lds.nphase+mod(r,lds.nphase)+1;
r = (0:((lds.ntst+1)*lds.nphase-1));
lds.multi_r2 = floor(r./lds.nphase)*lds.ncol_coord+mod(r,lds.nphase)+1;
lds.BranchParam = lds.ActiveParams;
% ------------------------------------------------------
function WorkspaceDone(x,v,s)


%------------------------------------------------------------------
function jac = BVP_BPC_jac(BVP_jac_f,BVP_jac_bc1,x,p,T,npar,nc)
global lds

p = num2cell(p);
pars1 = lds.ncoords+1;
pars2 = lds.ncoords+2;
pars3 = lds.ncoords+3;

jac = zeros(lds.ncoords+2,lds.ncoords+3);
ups = reshape(x,lds.nphase,lds.tps);
% function
range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
for j=lds.tsts
  xp = ups(:,range0)*lds.wt;
  jac(range1,[range2 pars1 pars2 pars3]) = feval(BVP_jac_f,lds.func,xp,p,T,j);
  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end
% boundary conditions
range  = (lds.tps-1)*lds.nphase+ (lds.phases);
range1 = lds.ncoords-lds.nphase+lds.phases;
jac(range,[lds.phases range1 pars1 pars2 pars3]) = feval(BVP_jac_bc1);

% integral constraint
jac(range(end)+1,[1:lds.ncoords pars3])= [BVP_LC_jac_ic lds.BPC_psi(end)];

range=range(end)+2:3;
jac(range,1:lds.ncoords+2) = lds.BPC_phi;

%------------------------------------------------------------

%SD:continues limit cycle of odefile
