function out = limitpointcycle
%
% LPCycles curve definition file
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
%--------------------------------------------------
function func = curve_func(arg)
global lds 
  [x,p,T] = rearr(arg);
  f = BVP('BVP_LC_f','BVP_LC_bc','BVP_LC_ic',x,p,T);
  % append g
  J = BVP_LPCjac2('BVP_LPC_f','BVP_LPC_bc1','BVP_LPC_bc2',x,p,T,1,1);
  
  b = []; b(lds.ncoords+2)=1;
  if lds.LPC_switch == 1
      s = J'\b';
      lds.LPC_new_psi = s(1:lds.ncoords+1)';
  else
      s = J\b';
      lds.LPC_new_phi = s(1:lds.ncoords+1)';
  end  
  G = s(lds.ncoords+2);
  f(end+1) = G;
  func = f;
  
%---------------------------------------------------------
function jac = jacobian(varargin)
global lds 
  [x,p,T] = rearr(varargin{1});
  jac = BVP_jac('BVP_LCX_jac',x,p,T,3,2);
  
  % calculate v and w
  fx = varargin{1};
  j = size(jac,1)+1;
  jac(j,:) = 0;
  b = []; b(lds.ncoords+2)=1;b=b';
  [x,p,T] = rearr(fx);
  J = BVP_LPCjac2('BVP_LPC_f','BVP_LPC_bc1','BVP_LPC_bc2',x,p,T,1,1);

  sn = J\b;
  st = J'\b;
  v = sn(1:lds.ncoords)';
  S = sn(lds.ncoords+1);
  w1 = st(1:lds.ncoords-lds.nphase)';
  w3 = st(end-1);

  % calculate g'
  ups = reshape(x,lds.nphase,lds.tps);
  p = num2cell(p);

  
  range0 = lds.cols_p1;
  range1 = lds.col_coords;
  range2 = lds.cols_p1_coords;
  
  t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
  kr1 = fix(t/lds.nphase);
  kr2 = rem(t,lds.nphase)+1;
  J = zeros(1,lds.ncoords+3);
  J(1,[lds.coords lds.ncoords+(2:3)]) = BVP_jac_LPC_ic(ups,p,v');

  for tstpt=lds.tsts
    xp  = ups(:,range0)*lds.wt;
    cv  = v(range2)';
    cw1 = w1(range1);
   
    range = lds.phases;
    for c=lds.cols
      xt = xp(:,c);
      sysj  = cjac(lds.func,lds.Jacobian,xt,p,lds.ActiveParams);
      sysp  = cjacp(lds.func,lds.JacobianP,xt,p,lds.ActiveParams);
      sysh  = chess(lds.func,lds.Jacobian,lds.Hessians,xt,p,lds.ActiveParams);
      syshp = chessp(lds.func,lds.Jacobian,lds.HessiansP,xt,p,lds.ActiveParams);
      wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
      for d=lds.phases
        sh(:,d) = (wtk.*sysh(:,kr2,d))*cv;
      end      
      t1 = T* wtk.*sh(:,kr2)+wtk.*sysj(:,kr2)*S;
      t2 = (wtk.*sysj(:,kr2))*cv;
      t3 = T* wtk.*syshp(:,kr2,1)* cv + sysp(:,1)*S;
      t4 = T* wtk.*syshp(:,kr2,2)* cv + sysp(:,2)*S;
      syshess(range,:) = [t1 t2 t3 t4];      
      range = range + lds.nphase;
    end
    jac(j,[range2 lds.ncoords+(1:3)]) = jac(j,[range2 lds.ncoords+(1:3)]) + cw1*syshess;
    
    range0 = range0 + lds.ncol;
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol_coord;
  end
  jac(j,[lds.coords lds.ncoords+(2:3)]) = jac(j,[lds.coords lds.ncoords+(2:3)])-w3*J(1,[lds.coords lds.ncoords+(2:3)]); 
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
    s.data.timemesh = lds.msh;
    s.data.ntst = lds.ntst;
    s.data.ncol = lds.ncol;
    s.data.parametervalues = p;
    s.data.T = T;
    s.data.multipliers=lds.multipliers;
    varargout{3} = s;
  end

  % update upoldp
  p1 = num2cell(p);
  for i=1:lds.tps
    lds.upoldp(:,i) = T*feval(lds.func, 0, lds.ups(:,i), p1{:});
  end
  
  % calculate multipliers if requested
  if lds.CalcMultipliers
    lds.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,varargin{1},[]));
    lds.multipliersX = varargin{1};
  end
  if lds.CalcMultipliers==0
      lds.multipliers = [];
  end
  varargout{2} = [lds.msh'; lds.multipliers];

  % all done succesfully
  varargout{1} = 0;

%---------------------------------------------------------
function option = options
global lds 
  % Check for symbolic derivatives in odefile
  
  symord = 0; 
  if ~isempty(lds.Jacobian), symord = 1; end
  if ~isempty(lds.Hessians), symord = 2; end
  if ~isempty(lds.Der3), symord = 3; end
  if ~isempty(lds.Der4), symord = 4; end
  if ~isempty(lds.Der5), symord = 5; end
  
  bp=size(lds.BranchParams,2);  
  option = contset;
  switch lds.nphase
      case {1 2}
          option=contset(option,'IgnoreSingularity',[bp+1 bp+2 bp+3 bp+4]);
      case 3
          option=contset(option,'IgnoreSingularity',[bp+3]);
  end
  option = contset(option, 'SymDerivative', symord);
  if lds.nphase < 4
      option=contset(option,'IgnoreSingularity',[bp+3 bp+4]);
  else
      [lds.bialt_M1,lds.bialt_M2,lds.bialt_M3,lds.bialt_M4]=bialtaa(lds.nphase-2);
  end
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', [0 0 0]);
  symordp = 0;
  if ~isempty(lds.JacobianP), symordp = 1; end
  if ~isempty(lds.HessiansP),  symordp = 2; end
  option = contset(option, 'SymDerivative', symord);
  option = contset(option, 'SymDerivativeP', symordp);
  
%------------------------------------------------------------------------
function [out, failed] = testf(id, x0, v)
global lds cds
 [x,p,T] = rearr(x0);
 %BP
bp = size(lds.BranchParams,2);     
if any(ismember(bp+3:bp+4,id)) & (isempty(lds.multipliersX)|(lds.multipliersX~=x0))
    lds.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,x0,[])); 
    lds.multipliersX = x0;
end

 failed = [];
   J = BVP_LPCjac2('BVP_LPC_f','BVP_LPC_bc1','BVP_LPC_bc2',x,p,T,1,1);

 b = []; b(lds.ncoords+2)=1;
 wext = J'\b';
 %out = zeros(1,lds.BranchParams+4);
 out=zeros(1,size(lds.BranchParams,2)+4);
 lastwarn('');
 if any(ismember(bp+1:bp+2,id))
     pt = num2cell(p);
     vext = J\b';
     vn = vext(1:lds.ncoords)';
     S = vext(lds.ncoords+1);
     v1 =(T/S)*vn;
     w1 = wext(1:lds.ncoords-lds.nphase)';
     
     ups = reshape( x,lds.nphase,lds.tps);
     vps = reshape(v1',lds.nphase,lds.tps);
     range1 = lds.cols_p1;
     range2 = lds.phases;
     range3 = lds.col_coords;
     range4 = lds.cols_p1_coords;
     t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
     kr1 = fix(t/lds.nphase);
     kr2 = rem(t,lds.nphase)+1;
     fxjac = spalloc(lds.ncoords-lds.nphase,lds.ncoords,(lds.ncol*lds.nphase)*((lds.ncol+1)*lds.nphase)*lds.ntst);
     fxhess = zeros(1,lds.ncoords-lds.nphase);
     v2 = zeros(1,lds.ncoords-lds.nphase);
     for j=lds.tsts
         % value of polynomial on each collocation point
         xp = ups(:,range1)*lds.wt;
         vp = vps(:,range1)*lds.wt;
         v3 = v1(range4);
         range5 = lds.phases;         
         sysjac = zeros(lds.ncol_coord,(lds.ncol+1)*lds.nphase);
         % evaluate function value on each collocation point
         for c=lds.cols                    
             xt = xp(:,c);
             wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
             hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);             
             sysjac(range5,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',cjac(lds.func,lds.Jacobian,xt,pt,lds.ActiveParams));
             v2(range2) = vp(:,c);
             for d=lds.phases
                 sh(:,d) = (wtk.*hess(:,kr2,d))*v3';
             end
             fxhess(range2)=wtk.*sh(:,kr2)*v3';
             range2 = range2+lds.nphase;
             range5 = range5+lds.nphase;
         end   
         fxjac(range3,range4)=sysjac;
         range1 = range1+lds.ncol;
         range3 = range3 + lds.ncol_coord;
         range4 = range4 + lds.ncol_coord;
     end
     res=w1*v2';
 end
 if ismember(bp+1,id)%R1
     out(bp+1) = res;
 end
 if ismember(bp+2,id)%CPC
     sysjacv1=fxjac*v1';
     out(bp+2)= (w1*(fxhess'+2*sysjacv1))/(2*res);
 end  
 if ismember(bp+3,id)%LPNS
    tmpmon = (lds.monodromy - eye(lds.nphase))^2;
    [Q,R,E] = qr(tmpmon);
    Q1 = Q(1:lds.nphase,1:lds.nphase-2);
    A = [Q1'*lds.monodromy*Q1 zeros(lds.nphase-2,1)];
    A = A(lds.bialt_M1).*A(lds.bialt_M2)-A(lds.bialt_M3).*A(lds.bialt_M4);
    A = A-eye(size(A,1));   
    out(bp+3) = det(A);
 end
 if ismember(bp+4,id)%LPPD
    A = lds.monodromy;
    A = A + eye(size(A,1));
    out(bp+4) = det(A);
 end

 if any(ismember(1:bp,id))%BPC
     wex = wext(1:lds.ncoords)';
     bp = size(lds.BranchParams,2);
     jac = BVP_jacbr(x,p,T,bp);
     id(id==(bp+1))=[];
     id(id==(bp+2))=[];
     id(id==(bp+3))=[];
     id(id==(bp+4))=[];
     out(id) = wex*jac(:,id);
 end
 if ~isempty(lastwarn)
     failed = [failed id];
 end            
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

%-------------------------------------------------------------
function [failed,s] = process (id, x, v, s)
global lds 
 bp = size(lds.BranchParams,2);
 
  % WM: Removed SL array
  switch id
      case {bp+1}
        s.data.c = nf_R1(x); 
        fprintf('Resonance 1:1 (period = %e, parameters = %e, %e)\n',x(length(x)-2),x(length(x)-1),x(length(x)));
        fprintf('ab=%d\n',s.data.c); 
        s.msg  = sprintf('Resonance 1:1');
      case {bp+2}
          s.data.cpccoefficient = nf_CPC(x);
          fprintf('label = CPC (period = %e, parameters = %e, %e)\n',x(length(x)-2),x(length(x)-1),x(length(x))); 
          fprintf('c = %d\n', s.data.cpccoefficient);          
          s.msg = sprintf('Cusp point'); 
      case {bp+3}
          s.data.c = nf_LPNS(x);                     
          fprintf('label = LPNS (period = %e, parameters = %e, %e)\n',x(length(x)-2),x(length(x)-1),x(length(x)));
          fprintf('(s,theta,E)=(%d, %d, %d)\n',s.data.c); 
          s.msg = sprintf('Fold-Neimark Sacker');
      case {bp+4}
          s.data.ffcoefficients = nf_FF(x);
          fprintf('Fold-flip (period = %e, parameters = %e, %e)\n',x(length(x)-2),x(length(x)-1),x(length(x)));      
          fprintf('(b11,a20,a02,cns) = (%d,%d,%d,%d)\n', s.data.ffcoefficients);
          s.msg = sprintf('Fold-flip');                    
      otherwise
          fprintf('label = %s ', s.label); 
          s.msg = sprintf('Branch point (parameter %d)',lds.BranchParams(id));
  end

  failed = 0;

%----------------------------------------------------------
function [S,L] = singmat
global lds 
  if isnan(lds.BranchParams)
    lds.BranchParams=[]; 
  end
  bp = size(lds.BranchParams,2);                    
  S = [];
  for i = 1:bp+1
      S(i,end+1)=0;
      S(i,1:end-1)=8;
      S(1:end-1,i)=8;
  end
  S(bp+2,end+1)=0;
  S(bp+2,1:end-2)=8;
  S(bp+2,end-1)=8;
  S(1:end-1,bp+2)=8;
   S(bp+3,end+1)=0;
  S(bp+3,1:end-2)=8;
  S(bp+3,end-1)=8;
  S(1:end-1,bp+3)=8;
  S(bp+3,bp+1)=1;
   S(bp+4,end+1)=0;
  S(bp+4,1:end-2)=8;
  S(bp+4,end-1)=8;
  S(1:end-1,bp+4)=8;

  
  L = 'R1';
  for i = 1:bp
      L(i,1:4)= strcat('BPC',num2str(lds.BranchParams(i)));
  end
  L(bp+1,1:4)='R1  ';
  L(bp+2,1:4)='CPC ';
  L(bp+3,1:4)='LPNS';
  L(bp+4,1:4)='LPPD';
 
%----------------------------------------------------------
function [x,v] = locate( varargin)
%----------------------------------------------------------
function varargout = init(varargin)
  WorkspaceInit(varargin{1:2});
  % all done succesfully
  varargout{1} = 0;
%----------------------------------------------------------
function done
  WorkspaceDone;

% -------------------------------------------------------
function [res,x,v] = adapt(x,v)
global lds
% calculate phi and psi for next point
if lds.LPC_switch == 0
  lds.LPC_phi = lds.LPC_new_phi/norm(lds.LPC_new_phi);
else
  lds.LPC_psi = lds.LPC_new_psi/norm(lds.LPC_new_psi);
end
lds.LPC_switch = 1-lds.LPC_switch;

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
function jac = BVP_LPCjac2(BVP_jac_f,BVP_jac_bc1,BVP_jac_bc2,x,p,T,npar,nc)
global lds

p = num2cell(p);
pars1 = lds.ncoords+1;
pars2 = lds.ncoords+1+(1:npar);

jac = spalloc(lds.ncoords+2,lds.ncoords+2,(lds.ncol*lds.nphase)*((lds.ncol+1)*lds.nphase+2)*lds.ntst+3*lds.ncoords);
ups = reshape(x,lds.nphase,lds.tps);
% function
range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
for j=lds.tsts
  xp = ups(:,range0)*lds.wt;
  jac(range1,[range2 pars1 pars2]) = feval(BVP_jac_f,lds.func,xp,p,T,j);
  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end
% boundary conditions
range  = (lds.tps-1)*lds.nphase+ (lds.phases);
range1 = lds.ncoords-lds.nphase+lds.phases;
jac(range,[lds.phases range1 pars1 pars2]) = feval(BVP_jac_bc1);

range = range(lds.nphase)+1;
jac(range,[lds.coords end]) = feval(BVP_jac_bc2,lds.func,ups,p);


% integral constraint
range=range+1;
jac(range,1:lds.ncoords+1) = lds.LPC_phi;

% -------------------------------------------------------------

function jac = BVP_jac(BVP_func,x,p,T,pars,nc)
global lds 
 
p2 = num2cell(p);
jac = feval(BVP_func,lds.func,x,p,T,pars,nc,lds,p2,lds.Jacobian,lds.ActiveParams,lds.JacobianP); 
jac = sparse(full(jac));


% ---------------------------------------------------------------
function jac = BVP_jacbr(x,p,T,npar)
global lds

ups = reshape(x,lds.nphase,lds.tps);
p = num2cell(p);
pars = 1:npar;
jac = zeros(lds.ncoords,npar);

range0 = lds.cols_p1;
range1 = lds.col_coords;

for j=lds.tsts
  xp = ups(:,range0)*lds.wt;
  jac(range1, pars) = BVP_LC_jacbr_f(lds.func,xp,p,T,j);
  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
end



% ---------------------------------------------------------------


function WorkspaceInit(x,v)
global cds lds


lds.LPC_new_phi = lds.LPC_phi;
lds.LPC_new_psi = lds.LPC_psi;
lds.LPC_switch = 0;

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



