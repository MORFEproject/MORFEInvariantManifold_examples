function [x0,v0]= init_BP_BP(odefile, x, p, ap, bp)
%
% [x0,v0] = init_BP_BP(odefile, x, p, ap, bp)
%
% Initializes a branch point continuation from a BP point
% 
%
global cds bpds 

% check input
 if size(ap,2)~=3
     errordlg('Three active parameter are needed for a Branchpoint bifurcation curve continuation');
 end
 if size(bp,2)~=1
     errordlg('One branch parameter is needed for a Branchpoint bifurcation curve continuation');
 end
% initialize bpds
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

bpds.odefile = odefile;
bpds.func = func_handles{2};
bpds.Jacobian  = func_handles{3};
bpds.JacobianP = func_handles{4};
bpds.Hessians  = func_handles{5};
bpds.HessiansP = func_handles{6};
bpds.Der3 = func_handles{7};
bpds.Der4 = func_handles{8};
bpds.Der5 = func_handles{9};

siz = size(func_handles,2);
if siz > 9
    j=1;
    for i=10:siz
        bpds.user{j}= func_handles{i};
        j=j+1;
    end
end
bpds.nphase = size(x,1);
bpds.ActiveParams = ap;
bpds.BranchParam =  bp;
bpds.P0 = p;
x0 = x ; 
x0(bpds.nphase+(1:3),:) = p(ap);
if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end
cds.curve = @branchpoint;
cds.ndim = length(x)+3;
[x,p] =rearr(x0); p = num2cell(p);
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_options = curvehandles{3};
cds.curve_jacobian = curvehandles{4};
cds.curve_hessians = curvehandles{5};
cds.options = feval(cds.curve_options);
cds.options = contset(cds.options,'Increment',1e-5);
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.symjac = 1;
cds.symhess = 0;
[Q,R,E]=qr([cjac(bpds.func,bpds.Jacobian,x,p,bpds.ActiveParams) cjacbr(bpds.func,bpds.JacobianP,x,p,bpds.ActiveParams,bpds.BranchParam)]);
bpds.borders.w=Q(:,end);
R(end,end)=0;
bpds.borders.v= E*[R(1:end-2,1:end-2)\-R(1:end-2,end-1:end);eye(2)];
% calculate eigenvalues
v0=[];
cds = rmfield(cds,'options');


% ---------------------------------------------------------------
function [x,p] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global bpds
p = bpds.P0;
p(bpds.ActiveParams) = x0((bpds.nphase+1):end);
x = x0(1:bpds.nphase);



