function [x0,v0]= init_EP_EP(odefile, x, p, ap)
%
%[x0,v0]= init_EP_EP(odefile, x, p, ap)
%
% Defines odefile
% Sets all parameters for an equilibrium continuation (p)
% and the active parameter (ap)
%
global cds eds
eds=[];
x0 = [x;p(ap)];
v0 = [];
eds.P0 = p;
eds.ActiveParams = ap;
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
if ~isfield(cds,'options') || ~isfield(cds,'options')
    cds.options = contset();
end
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.options.Increment = 1e-5;
cds.ndim = length(x)+1;
cds.symjac = 1;
cds.symhess = 1;


eds.odefile = odefile;
eds.func = func_handles{2};
eds.Jacobian  = func_handles{3};
eds.JacobianP = func_handles{4};
eds.Hessians  = func_handles{5};
eds.HessiansP = func_handles{6};
eds.Der3 = func_handles{7};
eds.Der4 = func_handles{8};
eds.Der5 = func_handles{9};
siz = size(func_handles,2);
if siz > 9
    j=1;
    for i=10:siz
        eds.user{j}= func_handles{i};
        j=j+1;
    end
end
eds.v = [];

eds.nphase = length(x);
n = eds.nphase;
a = reshape(1:(n^2),n,n);
[bia,bin,bip] = bialt(a);
if any(any(bip))
    [eds.BiAlt_M1_I,eds.BiAlt_M1_J,eds.BiAlt_M1_V] = find(bip);
else
    eds.BiAlt_M1_I=1;eds.BiAlt_M1_J=1;eds.BiAlt_M1_V=n^2+1;
end    
if any(any(bin))
    [eds.BiAlt_M2_I,eds.BiAlt_M2_J,eds.BiAlt_M2_V] = find(bin);
else
     eds.BiAlt_M2_I=1;eds.BiAlt_M2_J=1;eds.BiAlt_M2_V=n^2+1;
end
if any(any(bia))
    [eds.BiAlt_M3_I,eds.BiAlt_M3_J,eds.BiAlt_M3_V] = find(bia);
else
    eds.BiAlt_M3_I=1;eds.BiAlt_M3_J=1;eds.BiAlt_M3_V=n^2+1;
end

jac = cjac(eds.func,eds.Jacobian,x,num2cell(p),eds.ActiveParams);
jac = jac(:,1:cds.ndim-1);
jac(:,end+1) = 0;
A1 = sparse(eds.BiAlt_M1_I,eds.BiAlt_M1_J,jac(eds.BiAlt_M1_V));
A2 = sparse(eds.BiAlt_M2_I,eds.BiAlt_M2_J,jac(eds.BiAlt_M2_V));
A3 = sparse(eds.BiAlt_M3_I,eds.BiAlt_M3_J,jac(eds.BiAlt_M3_V));
A = A1-A2+A3;
[Q,R,E] = qr(full(A));
eds.bigW = Q(:,end);
eds.bigV = E(:,end);
eds.bigD = 0;

