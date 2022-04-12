function [x,v] = init_PD_LC(odefile, x, s, ntst, ncol, h)
%
% [x0,v0] = init_PD_LC(odefile, x, s, ntst, ncol)
%
% Initializes a limit cycle continuation of a double period cycle
% from a Period Doubling bifurcation detected during a previous
% limit cycle continuation.
%
global lds cds

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
else
    lds.user=[];
end


lds.nphase = floor((size(x,1)-2)/(s.data.ntst*s.data.ncol+1));
lds.P0 = s.data.parametervalues;
lds.T = 2*s.data.T;
set_ntst_ncol(s.data.ntst,s.data.ncol,s.data.timemesh);

% get x
x = x(:,s.index);

% generate double period solution
w = s.data.phi/norm(s.data.phi);

n_par = size(lds.ActiveParams,2);
if n_par==2
    x(end-n_par+1:end)=lds.P0(lds.ActiveParams);
else
    x(end-n_par+1) = lds.T;
    x(end) = lds.P0(lds.ActiveParams);
end

 if n_par==1
    x = [x(1:(lds.ncoords-lds.nphase))+h*w(1:(lds.ncoords-lds.nphase))';x(lds.coords)-h*w';lds.T;x(end)];
 else
     x = [x(1:(lds.ncoords-lds.nphase))+h*w(1:(lds.ncoords-lds.nphase))';x(lds.coords)-h*w';lds.P0(lds.ActiveParams)];     
 end
v = [w(1:(lds.ncoords-lds.nphase))';-w';0;0];
set_ntst_ncol(2*s.data.ntst,s.data.ncol,[lds.msh./2 (lds.msh(2:end)+1)./2])

% generate a new mesh and interpolate
[x,v]=new_mesh(x,v,ntst,ncol);
