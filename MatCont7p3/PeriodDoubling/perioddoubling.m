function out = perioddoubling
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
global lds
  [x,p,T] = rearr(arg);
  func = BVP('BVP_LC_f','BVP_LC_bc','BVP_LC_ic',x,p,T);

  % append g
  J = BVP_jac('BVP_PD_jac',x,p,T,1,1);
  b = []; b(lds.ncoords+1)=1;
  if lds.PD_switch == 1
    s = J'\b';
    lds.PD_new_psi = reshape(s(lds.coords),lds.nphase,lds.tps);
  else
    s = J\b';
    lds.PD_new_phi = reshape(s(lds.coords),lds.nphase,lds.tps);
  end
  G = s(lds.ncoords+1);
  func(end+1) = G;
  
%-------------------------------------------------------------
function jac = jacobian(varargin)
global lds 
%elseif strcmp(arg, 'jacobian')
  fx = varargin{1};
  [x,p,T] = rearr(fx);
  jac = BVP_jac('BVP_LCX_jac',x,p,T,3,2);
  % append g'
  
  % calculate v and w
  j = size(jac,1)+1;
  jac(j,:) = 0;
  b = []; b(lds.ncoords+1)=1;
  [x,p,T] = rearr(fx);
  J = BVP_jac('BVP_PD_jac',fx,p,T,1,1);
  sn = J\b';
  st = J'\b';
  v = st(lds.coords)';
  w = sn(lds.coords)';
  
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
    xp = ups(:,range0)*lds.wt;
    cw = w(range2)';
    cv = v(range1);
    
    range = lds.phases;
    for c=lds.cols
      xt = xp(:,c);
      sysj  = cjac(lds.func,lds.Jacobian,xt,p,lds.ActiveParams);
      sysh  = chess(lds.func,lds.Jacobian,lds.Hessians,xt,p,lds.ActiveParams);
      syshp = chessp(lds.func,lds.Jacobian,lds.HessiansP,xt,p,lds.ActiveParams);
      
      wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
      for d=lds.phases
        sh(:,d) = (wtk.*sysh(:,kr2,d))*cw;
      end
      syshess(range,:) = [(wtk.*sh(:,kr2)) (wtk.*sysj(:,kr2))*cw (wtk.*syshp(:,kr2,1))*cw (wtk.*syshp(:,kr2,2))*cw];
      
      range = range + lds.nphase;
    end
    jac(j,[range2 lds.ncoords+(1:3)]) = jac(j,[range2 lds.ncoords+(1:3)]) + cv*syshess;
    
    range0 = range0 + lds.ncol;
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol_coord;
  end
  jac(j,:) = T*jac(j,:);
  jac(j,lds.ncoords+1) = jac(j,lds.ncoords+1)/T;
  
%----------------------------------------------------------
function hessians(varargin)
%----------------------------------------------------------
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
    s.data=[];
    s.data.timemesh = lds.msh;
    s.data.ntst = lds.ntst;
    s.data.ncol = lds.ncol;
    s.data.parametervalues = p;
    s.data.T = T;
    s.data.multipliers=lds.multipliers;
    varargout{3} = s;
  end

  % update upoldp
  pt = num2cell(p);
  for i=1:lds.tps
    lds.upoldp(:,i) = T*feval(lds.func, 0, lds.ups(:,i), pt{:});
  end
  
  % calculate multipliers if requested
  if lds.CalcMultipliers %& (isempty(lds.multipliersX)|lds.multipliersX~=varargin{1})
    lds.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,varargin{1},[]));
    lds.multipliersX = varargin{1};
  end
  if lds.CalcMultipliers==0
      lds.multipliers = [];
  end
  varargout{2} = [lds.msh'; lds.multipliers];
  
  % all done succesfully
  varargout{1} = 0;
  
%----------------------------------------------------------
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
          option=contset(option,'IgnoreSingularity',[1 2 3 4]);
      case 2
          option=contset(option,'IgnoreSingularity',[1 2 4]);
      case 3
          option=contset(option,'IgnoreSingularity',[4]);
  end
  option = contset(option, 'SymDerivative', symord);
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', [0 0 0 0]);
  symordp = 0;
  if ~isempty(lds.JacobianP), symordp = 1; end
  if ~isempty(lds.HessiansP),  symordp = 2; end
  option = contset(option, 'SymDerivativeP', symordp);
  
%----------------------------------------------------------  
function [out, failed] = testf(id, x0, v)
global cds lds
[xt,p,T] = rearr(x0);
ups = reshape(xt,lds.nphase,lds.tps);
pt = num2cell(p);

%add border
failed=[];out(3)=0;
if ismember(4,id) & (isempty(lds.multipliersX) | (lds.multipliersX~=x0))
    lds.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,x0,[])); 
    lds.multipliersX = x0;
end

%J = BVP_jac2('BVP_PD_jac_f','BVP_PD_jac_bc','BVP_PD_jac_ic',xt,p,T,1,1);
J = BVP_jac('BVP_PD_jac',x0,p,T,1,1);
[LJ,UJ] = lu(J);
b = []; b(lds.ncoords+1)=1; b=b';
wext = LJ'\(UJ'\b);
vext = UJ\(LJ\b);
% %compute borders

v1 = vext(1:lds.ncoords);
vps = reshape(v1,lds.nphase,lds.tps);
% rescale vext
ficd = dot(vps,vps);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
v1 = v1/sqrt(ic);
% rescale wext
w1 = wext(1:lds.ncoords-lds.nphase)';
w2 = reshape(w1,lds.nphase,lds.tps-1);
ic2 = zeros(1,lds.ncoords);
range1 = lds.cols_p1_coords;
range2 = 1 : lds.ncol;
wt = lds.wt';
for j=lds.tsts
    pw =w2(:,range2)*wt;
    ic2(range1) = ic2(range1)+reshape(pw,1,(lds.nphase*(lds.ncol+1)));
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol;
end
ic3 = norm(w1,1);
w1=w1/ic3;
c1=1/ic3*ic2*v1;
if any(ismember([2 3],id))
    %calculate psi*
    % change boundary conditions
    range = (lds.tps-1)*lds.nphase + (lds.phases);
    J(range,[lds.phases range]) = BVP_LC1_jac_bc;
    % remove borders
    J(:,end)=[];J(end,:)=[];
    
    %compute borders
    %add borders
    J=[J lds.PD2_psi;lds.PD2_phi' 0];
    psi = J'\b;
    if lds.PD2_switch
        lds.PD2_new_psi=psi;
    else
        lds.PD2_new_phi=J\b;
    end
    lds.PD2_switch=1-lds.PD2_switch;
    
    ups = reshape(xt,lds.nphase,lds.tps);
    
    % function
    range1 = lds.cols_p1;
    range2 = lds.phases;
    range3 = lds.col_coords;
    range4 = lds.cols_p1_coords;
    t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
    kr1 = fix(t/lds.nphase);
    kr2 = rem(t,lds.nphase)+1;
    for j=lds.tsts
        % value of polynomial on each collocation point
        xp = ups(:,range1)*lds.wt;
        vp = vps(:,range1)*lds.wt;
        v3 = v1(range4);
        range5 = lds.phases;         
        % evaluate function value on each collocation point
        for c=lds.cols                    
            xt = xp(:,c);
            sten = zeros(lds.nphase,lds.nphase);
            wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
            sjac  = cjac(lds.func,lds.Jacobian,xt,pt,lds.ActiveParams);
            hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
            tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
            f1(range2) = feval(lds.func, 0,  xt, pt{:});
            sysjac(range5,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',sjac);
            
            v2(range2) = vp(:,c);
            for d1=lds.phases
                sh(:,d1) = (wtk.*hess(:,kr2,d1))*v3;
                for d2=lds.phases
                    stens(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v3;
                end    
                sten(:,d1)=sten(:,d1)+wtk.*stens(:,kr2,d1)*v3;
            end
            fxhess(range2)=wtk.*sh(:,kr2)*v3;
            fxtens(range2)=wtk.*sten(:,kr2)*v3;
            range2 = range2+lds.nphase;
            range5 = range5+lds.nphase;
        end   
        fxjac(range3,range4)=sysjac;
        range1 = range1+lds.ncol;
        range3 = range3 + lds.ncol_coord;
        range4 = range4 + lds.ncol_coord;
    end
    %rescale psi
    psi = psi(1:lds.ncoords-lds.nphase)';
    res2 = psi*f1';
    out(2) = res2;
end
    
lastwarn('');
if ismember(1,id)
    out(1)=c1;
end
if ismember(3,id)
    psi = psi*(1/(2*res2));           
    %computation of a
    a=(psi*fxhess');
    
    %computation of h2
    psi1 = reshape(psi,lds.nphase,lds.tps-1);
    ic = zeros(1,lds.ncoords);
    range1 = lds.cols_p1_coords;
    range2 = 1 : lds.ncol;
    wt = lds.wt';
    for j=lds.tsts
        pw =psi1(:,range2)*wt;
        ic(range1) = ic(range1)+reshape(pw,1,(lds.nphase*(lds.ncol+1)));
        range1 = range1 + lds.ncol_coord;
        range2 = range2 + lds.ncol;
    end
    %add integral constraint h2
    J(lds.ncoords+1,[lds.coords]) =ic;
    b = [];
    b(1:lds.ncoords-lds.nphase) = fxhess-2*a*f1; b(lds.ncoords+1) = 0;
    h2 = J\b';
    range1 = lds.cols_p1;
    range2 = lds.phases;
    range4 = lds.cols_p1_coords;
    for j=lds.tsts
        % value of polynomial on each collocation point
        xp = ups(:,range1)*lds.wt;
        v3 = v1(range4);
        h3 = h2(range4);
        range5 = lds.phases;         
        % evaluate function value on each collocation point
        for c=lds.cols                    
            xt = xp(:,c);    
            wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
            hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
            for d1=lds.phases
                sh(:,d1) = (wtk.*hess(:,kr2,d1))*v3;
            end
            fxhessh2(range2) = wtk.*sh(:,kr2)*h3;
            range2 = range2+lds.nphase;               
        end   
        range1 = range1+lds.ncol;
        range4 = range4 + lds.ncol_coord;
    end
    out(3) = c1*1/3*(w1*(1/T*fxtens'+3*fxhessh2'-6/T*a*fxjac*v1));
end
if ismember(4,id) %PDNS
    A = lds.monodromy;
    A = A(lds.bialt_M1).*A(lds.bialt_M2)-A(lds.bialt_M3).*A(lds.bialt_M4);
    A = A-eye(size(A,1));   
    out(4) = det(A);
end

if ~isempty(lastwarn)
    failed = [failed id];
end
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
function [failed,s] = process(id, x, v, s)  
global lds
switch id
  case 1
    s.data.c = nf_R2(x); 
    fprintf('Resonance 1:2 (period = %e, parameters = %e, %e)\n',x(length(x)-2),x(length(x)-1),x(length(x)));
    fprintf('(a,b)=(%d, %d)\n',s.data.c); 
    s.msg  = sprintf('Resonance 2:1'); 
  case 2
      s.data.ffcoefficients = nf_FF(x);
      fprintf('Fold-flip (period = %e, parameters = %e, %e)\n',x(length(x)-2),x(length(x)-1),x(length(x)));      
      fprintf('(b11,a20,a02,cns) = (%d,%d,%d,%d)\n', s.data.ffcoefficients);
      s.msg = sprintf('Fold-flip');
    case 3
      s.data.gpdcoefficient = nf_GPD(x);
      fprintf('Generalized period doubling(period = %e, parameters = %e, %e)\n',x(length(x)-2),x(length(x)-1),x(length(x)));
      fprintf('e = %d\n', s.data.gpdcoefficient);
      s.msg  = sprintf('generalized period doubling');           
  case 4
      s.data.c = nf_PDNS(x);       
      fprintf('Flip-Neimark Sacker (period = %e, parameters = %e, %e)\n',x(length(x)-2),x(length(x)-1),x(length(x)));
      fprintf('(p11,p22,theta,delta,sign(l1))=(%d, %d, %d, %d, %d)\n',s.data.c); 
      s.msg = sprintf('Flip-Neimark Sacker');
      s.data.multipliers = lds.multipliers;
end
failed = 0;
%--------------------------------------------------------
function [S,L] = singmat
  S = [ 0 8 8 8;
        8 0 8 8
        1 1 0 8
        1 1 8 0];
  L = ['R2  ';'LPPD';'GPD ';'PDNS'];

%--------------------------------------------------------
function locate(varargin)
%--------------------------------------------------------
function varargout = init(varargin)
  WorkspaceInit(varargin{1:2});
  % all done succesfully
  varargout{1} = 0;
%--------------------------------------------------------
function done
   WorkspaceDone;
%--------------------------------------------------------
function [res,x,v] = adapt(x,v)
global lds
% calculate phi and psi for next point
if lds.PD_switch == 0
  lds.PD_phi = lds.PD_new_phi;
  s = 1/sqrt(BVP_PD_jac_ic*lds.PD_phi(lds.coords)');
  t = 1/sqrt(s*norm(lds.PD_psi(lds.coords)));
else
  lds.PD_psi = lds.PD_new_psi;
  s = 1/norm(lds.PD_psi(lds.coords));
  t = sqrt(s*sqrt(BVP_PD_jac_ic*lds.PD_phi(lds.coords)'));
end
lds.PD_phi = (s/t)*lds.PD_phi;
lds.PD_psi = (s*t)*lds.PD_psi;
lds.PD_switch = 1-lds.PD_switch;

[x,v]=adapt_mesh(x,v);
res = 1;



%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%--------------------------------------------------------------
function jac = BVP_jac(BVP_func,x,p,T,pars,nc)
global lds 
 
p2 = num2cell(p);
jac = feval(BVP_func,lds.func,x,p,T,pars,nc,lds,p2,lds.Jacobian,lds.ActiveParams,lds.JacobianP); 
jac = sparse(full(jac));

% ---------------------------------------------------------------
function WorkspaceInit(x,v)
global cds lds

lds.pars = lds.ncoords+(1:3);
lds.tsts = 1:lds.ntst;
lds.cols = 1:lds.ncol;
lds.ntstcol = lds.ntst*lds.ncol;

lds.PD_new_phi = lds.PD_phi;
lds.PD_new_psi = lds.PD_psi;
[lds.PD2_phi,lds.PD2_psi]=initborders(x,v);
lds.PD_switch = 0;
lds.PD2_switch = 0;
lds.PD2_new_phi = lds.PD2_phi;
lds.PD2_new_psi = lds.PD2_psi;

[lds.bialt_M1,lds.bialt_M2,lds.bialt_M3,lds.bialt_M4]=bialtaa(lds.nphase);
lds.CalcMultipliers = contget(cds.options, 'Multipliers', 0);
lds.multipliersX = [];
lds.multipliers = nan;
lds.monodromy = [];

r = (0:(lds.ntst*lds.nphase-1));
lds.multi_r1 = (floor(r./lds.nphase)+1)*lds.ncol_coord-lds.nphase+mod(r,lds.nphase)+1;
r = (0:((lds.ntst+1)*lds.nphase-1));
lds.multi_r2 = floor(r./lds.nphase)*lds.ncol_coord+mod(r,lds.nphase)+1;

% ------------------------------------------------------

function WorkspaceDone

% -------------------------------------------------------
function [q,p] = initborders(x,v)
global lds 
%SD:continues period doubling bifurcation of odefile using minimal extended system
[xt,p,T] = rearr(x);
J = BVP_jac('BVP_PD_jac',x,p,T,1,1);
J(end,:)=rand(1,lds.ncoords+1);
J(:,end)=rand(lds.ncoords+1,1);
J(end,end)=0;
  b = []; b(lds.ncoords+1)=1;
  q=J\b'; q=q(1:end-1);q=q/norm(q);
  p=J'\b';p=p(1:end-1);p=p/norm(p);
