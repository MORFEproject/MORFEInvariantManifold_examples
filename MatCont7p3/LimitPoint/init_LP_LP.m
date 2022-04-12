function [x0,v0]= init_LP_LP(odefile, x, p, ap,varargin)
%
% [x0,v0] = init_LP_LP(odefile, x, p, ap) or
% [x0,v0] = init_LP_LP(odefile, x, p, ap, bp)
% Initializes a fold continuation from a LP point
% 
%
global cds lpds

% check input
if size(ap,2)~=2
    errordlg('Two active parameter are needed for a Limitpoint bifurcation curve continuation');
end
% initialize lpds
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


lpds.odefile = odefile;
lpds.func = func_handles{2};
lpds.Jacobian  = func_handles{3};
lpds.JacobianP = func_handles{4};
lpds.Hessians  = func_handles{5};
lpds.HessiansP = func_handles{6};
lpds.Der3 = func_handles{7};
lpds.Der4 = func_handles{8};
lpds.Der5 = func_handles{9};
siz = size(func_handles,2);
if siz > 9
    j=1;
    for i=10:siz
        lpds.user{j}= func_handles{i};
        j=j+1;
    end
end

lpds.nphase = size(x,1);
lpds.ActiveParams = ap;
lpds.P0 = p;
if size(varargin,1)>0,lpds.BranchParams=varargin{1};else lpds.BranchParams=[];end
x0 = x ; x0(lpds.nphase+1:lpds.nphase+2,:) = lpds.P0(ap);
cds.curve = @limitpoint;
cds.ndim = length(x)+2;
[x,p] =rearr(x0); p = num2cell(p);
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_options = curvehandles{3};
cds.curve_jacobian = curvehandles{4};
cds.curve_hessians = curvehandles{5};
cds.options = feval(cds.curve_options);
cds.symjac = 1;
cds.symhess = 0;

cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.options = contset(cds.options,'Increment',1e-5);
jac = cjac(lpds.func,lpds.Jacobian,x,p,lpds.ActiveParams);
% calculate eigenvalues
V = eig(jac);
[Y,i] = min(abs(V));
% ERROR OR WARNING
RED = jac-V(i)*eye(lpds.nphase);
lpds.borders.v = real(null(RED)); %HM: real because near BT the borders
lpds.borders.w = real(null(RED')); % might otherwise be complex.
v0=[];
%cds = rmfield(cds,'options');

%xxxxxxxxxxxxx
n = lpds.nphase;
a = reshape(1:(n^2),n,n);
[bia,bin,bip] = bialt(a);
if any(any(bip))
    [lpds.BiAlt_M1_I,lpds.BiAlt_M1_J,lpds.BiAlt_M1_V] = find(bip);
else
    lpds.BiAlt_M1_I=1;lpds.BiAlt_M1_J=1;lpds.BiAlt_M1_V=n^2+1;
end    
if any(any(bin))
    [lpds.BiAlt_M2_I,lpds.BiAlt_M2_J,lpds.BiAlt_M2_V] = find(bin);
else
     lpds.BiAlt_M2_I=1;lpds.BiAlt_M2_J=1;lpds.BiAlt_M2_V=n^2+1;
end
if any(any(bia))
    [lpds.BiAlt_M3_I,lpds.BiAlt_M3_J,lpds.BiAlt_M3_V] = find(bia);
else
    lpds.BiAlt_M3_I=1;lpds.BiAlt_M3_J=1;lpds.BiAlt_M3_V=n^2+1;
end

if lpds.nphase > 2
    jac(:,end+1) = 0;
    A1 = sparse(lpds.BiAlt_M1_I,lpds.BiAlt_M1_J,jac(lpds.BiAlt_M1_V));
    A2 = sparse(lpds.BiAlt_M2_I,lpds.BiAlt_M2_J,jac(lpds.BiAlt_M2_V));
    A3 = sparse(lpds.BiAlt_M3_I,lpds.BiAlt_M3_J,jac(lpds.BiAlt_M3_V));
    A = A1-A2+A3;
    [Q,R,E] = qr(full(A));
    lpds.bigW = Q(:,end);
    lpds.bigV = E(:,end);
    lpds.bigD = 0;
    lpds.adaptchoice = 1;
end


% ---------------------------------------------------------------
function [x,p] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global lpds

p = lpds.P0;
p(lpds.ActiveParams) = x0((lpds.nphase+1):end);
x = x0(1:lpds.nphase);



