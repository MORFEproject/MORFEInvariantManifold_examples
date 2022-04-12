function out = perioddoubling2
%
% Period doubling curve definition file
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
%---------------------------------------------------------
function func = curve_func(arg)

  [x,y,p,T] = rearr(arg);
  f = BVP('BVP_LC_f','BVP_LC_bc','BVP_LC_ic',x,p,T);

  % append second system
  f2 = BVP2('BVP_PD2_f','BVP_PD2_bc','BVP_PD2_ic',x,y,p,T);

  func = [f;f2];
  

%-------------------------------------------------------------
function jac = jacobian(varargin)

  [x,y,p,T] = rearr(varargin{1});
  f = BVP_jac('BVP_LC_jac_f','BVP_LC_jac_bc','BVP_LC_jac_ic',x,p,T,3,2);

  % append second systems derivatives  
  f2 = BVP2_jac('BVP_PD2_jac_f1','BVP_PD2_jac_f2','BVP_PD2_jac_bc','BVP_PD2_jac_ic',x,y,p,T,3,2);

  jac = [f; f2];

%----------------------------------------------------------
function hessians(varargin)
%----------------------------------------------------------
function varargout = defaultprocessor(varargin)
global lds cds
  [x,y,p,T] = rearr(varargin{1});
  [v,v2] = rearr(varargin{2});

  % update
  lds.ups = reshape(x,lds.nphase,lds.tps);
  lds.vps = reshape(v,lds.nphase,lds.tps);

  if nargin > 2
    % set data in special point structure
    s = varargin{3};
    s.data.timemesh = lds.msh;
    s.data.ntst = lds.ntst;
    s.data.ncol = lds.ncol;
    s.data.T = T;
    varargout{3} = s;
  end

  % update upoldp
  p = num2cell(p);
  for i=1:lds.tps
      lds.upoldp(:,i) = T*feval(lds.func, 0, lds.ups(:,i), p{:});
  end
  lds.vpoldp = reshape(y,lds.nphase,lds.tps);
  
  % calculate multipliers if requested
  if lds.CalcMultipliers && (lds.multipliersX~=varargin{1})
    lds.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,varargin{1},[]));
    lds.multipliersX = varargin{1};
  end
  
  % no special data
  varargout{2} = [lds.multipliers];

  % all done succesfully
  varargout{1} = 0;
    
%----------------------------------------------------------
function option = options
  % Check for symbolic derivatives in odefile
global lds 
  symord = 0; 
  if ~isempty(lds.Jacobian), symord = 1; end
  if ~isempty(lds.Hessians), symord = 2; end
  if ~isempty(lds.Der3), symord = 3; end
  if ~isempty(lds.Der4), symord = 4; end
  if ~isempty(lds.Der5), symord = 5; end

  option = contset;
  option = contset(option, 'SymDerivative', symord);
  option = contset(option, 'Singularities', 0);
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', [0 0 0]);
  symordp = 0;
  if ~isempty(lds.JacobianP), symordp = 1; end
  if ~isempty(lds.HessiansP),  symordp = 2; end
  option = contset(option, 'SymDerivativeP', symordp);

%----------------------------------------------------------  
function testf(varargin)
%----------------------------------------------------------
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
%--------------------------------------------------------
function process(varargin)  
%--------------------------------------------------------
function singmat(varargin)
%--------------------------------------------------------
function locate(varargin)
%--------------------------------------------------------
function varargout = init(varargin)
  WorkspaceInit(varargin{1:2});
  varargout{1} = 0;

%--------------------------------------------------------
function done
  WorkspaceDone;
%--------------------------------------------------------------
function [res,x,v] = adapt(x,v)
[x,v]=adapt_mesh2(x,v);
res = [];

% ---------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------


function [x,y,p,T] = rearr(x0)
global lds
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x), parameters (p) and period (T)
nap = length(lds.ActiveParams);

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+(1:nap));
x = x0(lds.coords);
y = x0(lds.coords+lds.ncoords);
T = x0(lds.PeriodIdx);

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

% -------------------------------------------------------------

function f = BVP2(BVP_f,BVP_bc,BVP_ic,x,y,p,T)
global lds

% extract ups
ups = reshape(x,lds.nphase,lds.tps);
vps = reshape(y,lds.nphase,lds.tps);
p = num2cell(p);

% function
range1 = lds.cols_p1;
range2 = lds.phases;
for j=lds.tsts
  xp = ups(:,range1)*lds.wt;
  yp = vps(:,range1)*lds.wt;
  t  = vps(:,range1)*lds.wpvec/lds.dt(j);
  for c=lds.cols
    f(range2) = feval(BVP_f,lds.func,t(:,c),xp(:,c),yp(:,c),p,T);
    range2 = range2+lds.nphase;
  end
  range1 = range1+lds.ncol;
end

% boundary conditions
f(range2) = feval(BVP_bc,vps(:,1),vps(:,lds.tps));
% integral constraint
f(lds.ncoords+1) = feval(BVP_ic,vps);

f = f';

% -------------------------------------------------------------

function jac = BVP_jac(BVP_jac_f,BVP_jac_bc,BVP_jac_ic,x,p,T,npar,nc)
global lds

ups = reshape(x,lds.nphase,lds.tps);
p = num2cell(p);
pars = 2*lds.ncoords+(1:npar);

jac = spalloc(lds.ncoords+1,lds.ncoords+1,(lds.ncol+4)*lds.nphase);
%jac = zeros(lds.ncoords+1,lds.ncoords+length(p)-1);

range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;

for j=lds.tsts
  xp = ups(:,range0)*lds.wt;

  jac(range1,[range2 pars]) = feval(BVP_jac_f,lds.func,xp,p,T,j);

  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end

% integral constraint
jac(lds.ncoords+1,[lds.coords]) = feval(BVP_jac_ic);

% boundary conditions
range = (lds.tps-1)*lds.nphase + (lds.phases);
jac(range,[lds.phases range lds.ncoords+(1:nc)]) = feval(BVP_jac_bc);
jac = sparse(full(jac));

% ---------------------------------------------------------------

function jac = BVP2_jac(BVP_jac_f1,BVP_jac_f2,BVP_jac_bc,BVP_jac_ic,x,y,p,T,npar,nc)
global lds

ups = reshape(x,lds.nphase,lds.tps);
vps = reshape(y,lds.nphase,lds.tps);
p = num2cell(p);
pars = 2*lds.ncoords+(1:npar);

jac = spalloc(lds.ncoords+1,lds.ncoords+1,(lds.ncol+4)*lds.nphase);
%jac = zeros(lds.ncoords+1,lds.ncoords+length(p)-1);

range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
range3 = lds.cols_p1_coords+(lds.ncoords);

for j=lds.tsts
  xp = ups(:,range0)*lds.wt;
  yp = vps(:,range0)*lds.wt;

  jac(range1,[range2 pars]) = feval(BVP_jac_f1,lds.func,xp,yp,p,T,j);
  jac(range1,range3) = feval(BVP_jac_f2,lds.func,xp,yp,p,T,j);

  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
  range3 = range3 + lds.ncol_coord;
end

% integral constraint
jac(lds.ncoords+1,lds.coords+lds.ncoords) = feval(BVP_jac_ic);

% boundary conditions
range = (lds.tps-1)*lds.nphase + (lds.phases);
jac(range,[lds.phases range lds.ncoords+(1:nc)]+lds.ncoords) = feval(BVP_jac_bc);


% ---------------------------------------------------------------


function WorkspaceInit(x,v)
global cds lds

lds.cols_p1 = 1:(lds.ncol+1);
lds.cols_p1_coords = 1:(lds.ncol+1)*lds.nphase;
lds.ncol_coord = lds.ncol*lds.nphase;
lds.col_coords = 1:lds.ncol*lds.nphase;
lds.coords = 1:lds.ncoords;
lds.pars = 2*lds.ncoords+(1:3);
lds.tsts = 1:lds.ntst;
lds.cols = 1:lds.ncol;
lds.phases = 1:lds.nphase;
lds.ntstcol = lds.ntst*lds.ncol;
lds.PeriodIdx = 2*lds.ncoords+1;

lds.idxmat = reshape(fix((1:((lds.ncol+1)*lds.ntst))/(1+1/lds.ncol))+1,lds.ncol+1,lds.ntst);
lds.dt = lds.msh(lds.tsts+1)-lds.msh(lds.tsts);

lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);

lds.CalcMultipliers = contget(cds.options, 'Multipliers', 0);
lds.multipliersX = [];
lds.multipliers = nan;

r = (0:(lds.ntst*lds.nphase-1));
lds.multi_r1 = (floor(r./lds.nphase)+1)*lds.ncol_coord-lds.nphase+mod(r,lds.nphase)+1;
r = (0:((lds.ntst+1)*lds.nphase-1));
lds.multi_r2 = floor(r./lds.nphase)*lds.ncol_coord+mod(r,lds.nphase)+1;

% ------------------------------------------------------

function WorkspaceDone
%SD:continues period doubling bifurcation of odefile using full extended system
