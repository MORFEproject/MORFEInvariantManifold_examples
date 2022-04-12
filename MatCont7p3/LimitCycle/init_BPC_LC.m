function [x,v] = init_BPC_LC(odefile, x, v, s, ntst, ncol, h)
%
% [x0,v0] = init_BPC_LC(odefile, x, v, s, ap, ntst, ncol, h)
%
% Initializes a secondary limit cycle continuation from a branch point of cycles calculated
% in a previous run. (it does branch switching)
%
global lds cds
% check input
% initialize lds
func_handles = feval(odefile);
symord = 0; 
symordp = 0;

if     ~isempty(func_handles{9}),   symord = 5; 
elseif ~isempty(func_handles{8}),   symord = 4; 
elseif ~isempty(func_handles{7}),   symord = 3; 
elseif ~isempty(func_handles{5}),   symord = 2; 
elseif ~isempty(func_handles{3}),   symord = 1; 
end
if     ~isempty(func_handles{6}),   symordp = 2; 
elseif ~isempty(func_handles{4}),   symordp = 1; 
end
if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.symjac = 1;
cds.symhess = 0;


lds.odefile = odefile;
lds.func = func_handles{2};
lds.Jacobian  = func_handles{3};
lds.JacobianP = func_handles{4};
lds.Hessians  = func_handles{5};
lds.HessiansP = func_handles{6};
lds.Der3 = func_handles{7};
lds.Der4 = func_handles{8};
lds.Der5 = func_handles{9};
siz = size(func_handles,2);
if siz > 9
    j=1;
    for i=10:siz
        lds.user{j}= func_handles{i};
        j=j+1;
    end
end

lds.nphase = (size(x,1)-2)/(s.data.ntst*s.data.ncol+1);
lds.BranchParam = lds.ActiveParams;
lds.ActiveParams = lds.ActiveParams;
lds.P0 = s.data.parametervalues;
set_ntst_ncol(s.data.ntst,s.data.ncol,s.data.timemesh);
% get x and v
n_par = size(lds.ActiveParams,2);
x0 = x(:,s.index);
oldv = v(:,s.index);
if n_par==2
    x0(end-n_par+1:end)=lds.P0(lds.ActiveParams);
else
    x0(end-n_par+1) = lds.T;
    x0(end) = lds.P0(lds.ActiveParams);
end
% generate a new mesh and interpolate

lds.cols_p1 = 1:(lds.ncol+1);
lds.cols_p1_coords = 1:(lds.ncol+1)*lds.nphase;
lds.ncol_coord = lds.ncol*lds.nphase;
lds.col_coords = 1:lds.ncol*lds.nphase;

lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);
lds.phases = 1:lds.nphase;

[xt,p,T] = rearr(x0);
ups = reshape(xt,lds.nphase,lds.tps);
p = num2cell(p);
pars1 = lds.ncoords+1;
pars2 = lds.ncoords+2;
jac = spalloc(lds.ncoords+1,lds.ncoords+2,(lds.ncol+4)*lds.nphase);

% function
range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
for j=lds.tsts
  xp = ups(:,range0)*lds.wt;
  jac(range1,[range2 pars1 pars2]) = bordBVP_BPC_f(lds.odefile,xp,p,T,j);
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
  p = lds.dt(j)*(ups(:,range1).*lds.pwi);
  ic(range2) = ic(range2)+p(lds.cols_p1_coords);
  range1 = range1 + lds.ncol;
  range2 = range2 + lds.ncol_coord;
end

jac(range(end)+1,1:lds.ncoords)= ic;
%compute borders
[Q,R,E] = qr(full(jac));
R(end,end) = 0;R(end,end-1) = 0;
pp = E*[R(1:end-1,1:end-2)\-R(1:end-1,end-1:end);eye(2)];
clear Q R E;
x1 = pp(:,1); x2 = pp(:,2);
alpha = -oldv'*x2;
beta  = oldv'*x1;
vnew  = alpha*x1+beta*x2;
v0 = vnew/norm(vnew);
[x,v] = new_mesh(x0,v0,s.data.ntst,s.data.ncol);
x = x + h * v;


% ---------------------------------------------------------------
function [x,p,T] = rearr(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x), parameters (p) and period (T)
global lds

nap = length(lds.ActiveParams);

p = lds.P0;
p(lds.ActiveParams) = x0(lds.PeriodIdx+(1:nap));
x = x0(lds.coords);
T = x0(lds.PeriodIdx);
