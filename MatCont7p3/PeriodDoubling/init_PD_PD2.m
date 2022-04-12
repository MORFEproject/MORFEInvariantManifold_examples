function [x,v] = pdinit_cycle2(odefile, x, s, ap, ntst, ncol)
%
% [x0,v0] = pdinit_cycle(odefile, x, v, s, ap, ntst, ncol)
%
% check input
if size(ap)~=[1 2]
    error('Two active parameters are needed for a Period Doubling bifurcation curve continuation');
end
global cds lds

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
lds.ActiveParams = ap;
lds.P0 = s.data.parametervalues;
set_ntst_ncol(s.data.ntst,s.data.ncol,s.data.timemesh);

% get x
x = x(:,s.index);

% append extra parameter
x(lds.ncoords+(2:3)) = lds.P0(ap);

% generate a new mesh and interpolate
hulp=s.data.phi';
[x,hulp]=new_mesh(x,hulp,ntst,ncol);

% append variational part to x
vps = reshape(hulp,lds.nphase,lds.tps);
lds.vpoldp = vps;
w = 1/sqrt(BVP_PD2_ic(vps)+1);
lds.vpoldp = w*vps;
x = [x(lds.coords);w*hulp;x(lds.ncoords+(1:3))];

% let continuer calculate tangent vector
v = [];