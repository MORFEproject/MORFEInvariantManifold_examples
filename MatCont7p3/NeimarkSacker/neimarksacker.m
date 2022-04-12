function out = neimarksacker
%
% NS curve definition file
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
%----------------------------------------------------
function func = curve_func(arg)
global lds 
  [x,p,T,k] = rearr(arg);
  
  f = BVP('BVP_LC_f','BVP_LC_bc','BVP_LC_ic',x,p,T);
  % append g
  J = BVP_NSjac('BVP_NS_jac',x,p,T,k,1,1);
  b = zeros((2*lds.tps-1)*lds.nphase+2,2);
  b(end-1:end,:) = eye(2);  
  s = J\b;  
  f(end+1) = s(end+lds.index1(1)-2,lds.index1(2));
  f(end+1) = s(end+lds.index2(1)-2,lds.index2(2));
  func = f;

%--------------------------------------------------------
function jac = jacobian(varargin)
global lds 
  [x,p,T,k] = rearr(varargin{1});

  jac = BVP_jac('BVP_LCX_jac',x,p,T,3,2);
  %append g'

  % calculate v and w

  j = size(jac,1)+1;
  jac(end,end+1) = 0;
  jac(j,:) = 0; jac(j+1,:) = 0;
  
  b = zeros((2*lds.tps-1)*lds.nphase+2,2);
  b(end-1:end,:) = eye(2);
  J = BVP_NSjac('BVP_NS_jac',x,p,T,k,1,1);
  sn = J\b;
  st = J'\b;
  v = sn(1:end-2,:)';
  w = st(1:end-2,:)';
  w1 = w(:,end-lds.nphase+1:end);
  
  % calculate g'
  ups = reshape(x,lds.nphase,lds.tps);
  p = num2cell(p);

  
  range1 = lds.col_coords;
  range2 = lds.cols_p1_coords;
  cv=[];
  
  t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
  kr1 = fix(t/lds.nphase);
  kr2 = rem(t,lds.nphase)+1;
  for ns = 1:2
      range3 = lds.cols_p1_coords;
      v1     = cv  ;
      range0 = lds.cols_p1;
      for tstpt = lds.tsts
          xp  = ups(:,range0)*lds.wt;
          cv = v(:,range2)';
          cw = w(:,range1);
          range = lds.phases;
          for c = lds.cols
              xt    = xp(:,c);
              sysj  = cjac(lds.func,lds.Jacobian,xt,p,lds.ActiveParams);
              sysh  = chess(lds.func,lds.Jacobian,lds.Hessians,xt,p,lds.ActiveParams);
              syshp = chessp(lds.func,lds.Jacobian,lds.HessiansP,xt,p,lds.ActiveParams);
              wtk   = lds.wt(kr1,c(ones(1,lds.nphase)))';
              for d = lds.phases
                  sh1(:,d) = (wtk.*sysh(:,kr2,d))*cv(:,lds.index1(2));
                  sh2(:,d) = (wtk.*sysh(:,kr2,d))*cv(:,lds.index2(2));
              end      
              t11 = T* wtk.*sh1(:,kr2);
              t21 = T* wtk.*sh2(:,kr2);
              t12 = (wtk.*sysj(:,kr2))*cv(:,lds.index1(2));
              t22 = (wtk.*sysj(:,kr2))*cv(:,lds.index2(2));
              t13 = T* wtk.*syshp(:,kr2,1)* cv(:,lds.index1(2));
              t23 = T* wtk.*syshp(:,kr2,1)* cv(:,lds.index2(2));
              t14 = T* wtk.*syshp(:,kr2,2)* cv(:,lds.index1(2));
              t24 = T* wtk.*syshp(:,kr2,2)* cv(:,lds.index2(2));
              syshess1(range,:) = [t11 t12 t13 t14];      
              syshess2(range,:) = [t21 t22 t23 t24];      
              range = range + lds.nphase;    
          end
          jac(j,[range3 lds.ncoords+(1:3)])   = jac(j,[range3 lds.ncoords+(1:3)])   + cw(lds.index1(1),:)*syshess1;
          jac(j+1,[range3 lds.ncoords+(1:3)]) = jac(j+1,[range3 lds.ncoords+(1:3)]) + cw(lds.index2(1),:)*syshess2;
          range0 = range0 + lds.ncol;
          range1 = range1 + lds.ncol_coord;
          range2 = range2 + lds.ncol_coord;
          range3 = range3 + lds.ncol_coord;
      end
  end
  jac(j,lds.ncoords+4)   = 2*w1(lds.index1(1),:)*v1(end-lds.nphase+1:end,lds.index1(2));
  jac(j+1,lds.ncoords+4) = 2*w1(lds.index2(1),:)*v1(end-lds.nphase+1:end,lds.index2(2));
%------------------------------------------------------------
function hessians(varargin)
%------------------------------------------------------------
function varargout = defaultprocessor(varargin)
global lds cds
  [x,p,T,k] = rearr(varargin{1});
  v = rearr(varargin{2});

  % update
  lds.ups = reshape(x,lds.nphase,lds.tps);
  lds.vps = reshape(v,lds.nphase,lds.tps);

   % update upoldp
  p1 = num2cell(p);
  for i=1:lds.tps
    lds.upoldp(:,i) = T*feval(lds.func, 0, lds.ups(:,i),  p1{:});
  end
  % calculate multipliers if requested
  if lds.CalcMultipliers %& (isempty(lds.multipliersX)|lds.multipliersX~=varargin{1})
    lds.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,varargin{1},[]));
    lds.multipliersX = varargin{1};
  end
  
  if nargin > 2
    % set data in special point structure
    s = varargin{3};
    if ~lds.CalcMultipliers | isempty(lds.multipliersX) | (lds.multipliersX~=varargin{1})
        s.data.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,varargin{1},[]));
        lds.multipliersX = varargin{1};
    else
        s.data.multipliers = lds.multipliers;
    end
    s.data.timemesh = lds.msh;
    s.data.ntst = lds.ntst;
    s.data.ncol = lds.ncol;
    s.data.parametervalues = p;
    s.data.T = T;
    s.data.multipliers=lds.multipliers;
    varargout{3} = s;
  end

  if lds.CalcMultipliers==0
      lds.multipliers = [];
  end
  varargout{2} = [lds.msh'; lds.multipliers];
  
  % all done succesfully
  varargout{1} = 0;
%------------------------------------------------------------    
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
  option = contset(option, 'SymDerivative', symord);
  %options = contset(options, 'Singularities', 0);
  if lds.nphase < 4
      option=contset(option,'IgnoreSingularity',[7 8]);
  else
      [lds.bialt_M1,lds.bialt_M2,lds.bialt_M3,lds.bialt_M4]=bialtaa(lds.nphase-2);
  end
  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', [0 0 0]);
  symordp = 0;
  if ~isempty(lds.JacobianP), symordp = 1; end
  if ~isempty(lds.HessiansP),  symordp = 2; end
  option = contset(option, 'SymDerivativeP', symordp);
  
%---------------------------------------------------------
function [out, failed] = testf(id, x0, v)
global lds cds

[xt,p,T,k] = rearr(x0);
failed = []; out(8)=0;
if any(ismember([7 8],id)) & (isempty(lds.multipliersX) | (lds.multipliersX~=x0))
    lds.multipliers = multipliers(cjac(cds.curve_func,cds.curve_jacobian,x0,[])); 
    lds.multipliersX = x0;
end
if any(ismember([5 6],id))
    pt = num2cell(p);
    theta = acos(k);
    %computing vext and wext
    [Jjac,Jjac2] = sBVP_jac('bordBVP_PD_jac_f','BVP_LC1_jac_bc',xt,p,T,theta,1,1);
    %add borders
    Jjac = [Jjac lds.NS1_psi;lds.NS1_phi.' 0];
    [LJ,UJ] = lu(Jjac);
    b = []; b(lds.ncoords+1)=1; b=b';
    vext = UJ\(LJ\b);
    wext = LJ'\(UJ'\b);
    if lds.NS1_switch
        lds.NS1_new_psi=wext;
    else
        lds.NS1_new_phi=vext;
    end
    lds.NS1_switch=1-lds.NS1_switch;
    
    %rescale vext
    
    v1 = vext(1:lds.ncoords);
    vps = reshape(v1,lds.nphase,lds.tps);
    ups = reshape(xt,lds.nphase,lds.tps);
    ficd = dot(vps,vps);
    ficdmat = ficd(lds.idxmat);
    ic = sum(lds.dt.*(lds.wi*ficdmat));
    v1 = v1/sqrt(ic);
    
    % rescale wext
    w1 = wext(1:lds.ncoords-lds.nphase).';
    ic3 = norm(w1,1);
    w1 = w1/ic3;
    
    %calculate psi*
    
    %add borders
    Jjac2=[Jjac2 lds.NS2_psi;lds.NS2_phi' 0];
    psi = Jjac2'\b;
    if lds.NS2_switch
        lds.NS2_new_psi=psi;
    else
        lds.NS2_new_phi=Jjac2\b;
    end
    lds.NS2_switch=1-lds.NS2_switch; 
    
    % function
    range1 = lds.cols_p1;
    range2 = lds.phases;
    range3 = lds.col_coords;
    range4 = lds.cols_p1_coords;
    range6 = lds.cols;
    t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
    kr1 = fix(t/lds.nphase);
    kr2 = rem(t,lds.nphase)+1;
    vps = reshape(v1,lds.nphase,lds.tps);
    
    for jk = lds.tsts
        % value of polynomial on each collocation point
        xp(:,range6) = ups(:,range1)*lds.wt;
        vp(range3) = vps(:,range1)*lds.wt;
        v3 = v1(range4);
        v4 = conj(v3);
        range5 = lds.phases;         
        % evaluate function value on each collocation point
        for c=lds.cols                    
            xt = xp(:,range6(c));
            sten = zeros(lds.nphase,lds.nphase);
            wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
            sjac  = cjac(lds.func,lds.Jacobian,xt,pt,lds.ActiveParams);
            hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
            tens = ctens3(lds.func,lds.Jacobian,lds.Hessians,lds.Der3,xt,pt,lds.ActiveParams);
            f1(range2) = feval(lds.func, 0,  xt, pt{:});
            sysjac(range5,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',sjac);        
            for d1=lds.phases
                sh(:,d1) = (wtk.*hess(:,kr2,d1))*v3;
                for d2=lds.phases
                    stens(:,d2,d1) = (wtk.*tens(:,kr2,d1,d2))*v3;
                end    
                sten(:,d1)=sten(:,d1)+wtk.*stens(:,kr2,d1)*v3;
            end
            fxhess(range2) = wtk.*sh(:,kr2)*v3;
            fxhessc(range2)= wtk.*sh(:,kr2)*v4;
            fxtens(range2) = wtk.*sten(:,kr2)*v4;
            range2 = range2 + lds.nphase;
            range5 = range5 + lds.nphase;
        end   
        fxjac(range3,range4) = sysjac;
        range1 = range1 + lds.ncol;
        range3 = range3 + lds.ncol_coord;
        range4 = range4 + lds.ncol_coord;
        range6 = range6 + lds.ncol;
    end
    alpha = dot(w1,vp);
    %rescale psi
    psi = psi(1:lds.ncoords-lds.nphase)';
    hres = psi*f1';
end
%lastwarn('');
if ismember(1,id)%  case 1 % R1
    out(1) = k-1;
end
if ismember(2,id)% case 2 % R2
    out(2) = k+1;
end
if ismember(3,id)%  case 3 % R3
    out(3) = k+1/2;
end
if ismember(4,id)%  case 4 % R4
    out(4) = k;
end
if ismember(5,id)%  case 5  %LPNS
    out(5)= hres;
end 
if ismember(6,id)% %case 6 %CH
    if abs(k) < 1
        psi = psi/hres;
        %computation of a
        a = dot(psi,real(fxhessc));    
        %calculate h20

        range0 = lds.cols_p1;
        range1 = lds.col_coords;
        range2 = lds.cols_p1_coords;
        range3 = lds.cols;
        twitheta = 2*1i*theta*lds.pwwt;
        for jk  = lds.tsts
            % evaluate part of Jacobian matrix
            Jjac(range1,range2) = bordBVP_PD_jac_f(lds.func,xp(:,range3),pt,T,jk)+twitheta;  
            range0 = range0 + lds.ncol;
            range1 = range1 + lds.ncol_coord;

            range2 = range2 + lds.ncol_coord;
            range3 = range3 + lds.ncol;
        end
        % remove borders
        Jjac(:,end)=[];Jjac(end,:)=[];
        b=[];
        b(1:lds.ncoords-lds.nphase) = fxhess; b(lds.ncoords) = 0;
        h20 =Jjac\b.';
        %computation of h11
        psi1 = reshape(psi,lds.nphase,lds.tps-1);
        ic = zeros(1,lds.ncoords);
        range1 = lds.cols_p1_coords;
        range2 = lds.cols;
        wt = lds.wt';
        for jk=lds.tsts        
            pw = psi1(:,range2)*wt;
            ic(range1) = ic(range1)+reshape(pw,1,(lds.nphase*(lds.ncol+1)));
            range1 = range1 + lds.ncol_coord;
            range2 = range2 + lds.ncol;
        end
        %add integral constraint h11
        Jjac2(lds.ncoords+1,[lds.coords]) =ic;
        b = []; 
        b(1:lds.ncoords-lds.nphase) = real(fxhessc)-a*f1; b(lds.ncoords+1) = 0;
        h11 = Jjac2\b.';
        range1 = lds.cols_p1;
        range2 = lds.phases;
        range4 = lds.cols_p1_coords;
        range6 = lds.cols;
        for jk=lds.tsts

            % value of polynomial on each collocation point
            v2 = v1(range4);
            v3 = conj(v1(range4));
            h113 = h11(range4);
            h203 = h20(range4);

            % evaluate function value on each collocation point
            for c=lds.cols                    
                xt = xp(:,range6(c));
                wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
                hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,pt,lds.ActiveParams);
                for d1=lds.phases
                    sh11(:,d1) = (wtk.*hess(:,kr2,d1))*h113;
                    sh20(:,d1) = (wtk.*hess(:,kr2,d1))*h203;
                end
                fxhessh11(range2) = wtk.*sh11(:,kr2)*v2;
                fxhessh20(range2) = wtk.*sh20(:,kr2)*v3;
                range2 = range2+lds.nphase;
            end   
            range1 = range1 + lds.ncol;
            range4 = range4 + lds.ncol_coord;
            range6 = range6 + lds.ncol;
        end
        out(6) = real(conj(alpha)*(1/2*dot(w1,(1/T)*(fxtens.')+2*(fxhessh11.')+(fxhessh20.'))-(a/T)*dot(w1,fxjac*v1))+abs(alpha)^2*(1i*a*theta)/(T^2));              
    else
        out(6) = 500;
    end
end
if ismember(7,id)%PDNS
    A = lds.monodromy;
    A = A + eye(size(A,1));
    out(7) = det(A);
end
if ismember(8,id)%NSNS
    tmpmon = lds.monodromy*lds.monodromy - 2*k*lds.monodromy + eye(lds.nphase);
    [Q,R,E] = qr(tmpmon);
    Q1 = Q(1:lds.nphase,1:lds.nphase-2);
    A = [Q1'*lds.monodromy*Q1 zeros(lds.nphase-2,1)];
    A = A(lds.bialt_M1).*A(lds.bialt_M2)-A(lds.bialt_M3).*A(lds.bialt_M4);
    A = A-eye(size(A,1));
    out(8) = det(A);
    lastwarn('')
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
function [failed,s] = process(id,x,v,s)
switch id
  case 1
    s.data.c = nf_R1(x); 
    fprintf('Resonance 1:1 (period = %e, parameters = %e, %e)\n',x(length(x)-3),x(length(x)-2),x(length(x)-1));
    fprintf('ab=%d\n',s.data.c);     
    s.msg  = sprintf('Resonance 1:1');
  case 2
    s.data.c = nf_R2(x); 
    fprintf('Resonance 1:2 (period = %e, parameters = %e, %e)\n',x(length(x)-3),x(length(x)-2),x(length(x)-1));
    fprintf('(a,b)=(%d, %d)\n',s.data.c); 
    s.msg  = sprintf('Resonance 1:2'); 
  case 3
    s.data.c = nf_R3(x); 
    fprintf('Resonance 1:3 (period = %e, parameters = %e, %e)\n',x(length(x)-3),x(length(x)-2),x(length(x)-1));
    fprintf('(b,Re(c))=(%d + (%d) i,%d)\n', [real(s.data.c(1)), imag(s.data.c(1)),s.data.c(2)]);
    s.msg  = sprintf('Resonance 1:3'); 
    
  case 4
    s.data.c = nf_R4(x); 
    fprintf('Resonance 1:4 (period = %e, parameters = %e, %e)\n',x(length(x)-3),x(length(x)-2),x(length(x)-1));
    fprintf('(A,d)=(%d+ (%d) i, %d+ (%d) i)\n',[real(s.data.c(1)), imag(s.data.c(1)),real(s.data.c(2)), imag(s.data.c(2))]); 
    s.msg  = sprintf('Resonance 1:4');    
  case 5
    s.data.c = nf_LPNS(x); 
    fprintf('Fold-Neimark Sacker bifurcation (period = %e, parameters = %e, %e)\n',x(length(x)-3),x(length(x)-2),x(length(x)-1));
    fprintf('(s,theta,E)=(%d, %d, %d)\n',s.data.c); 
    s.msg= sprintf('Fold-Neimark Sacker bifurcation');
  case 6
    s.data.c = nf_CH(x); 
    fprintf('Chenciner bifurcation (period = %e, parameters = %e, %e)\n',x(length(x)-3),x(length(x)-2),x(length(x)-1));
    fprintf('Re(e)=%d',s.data.c); 
    s.msg = sprintf('Chenciner bifurcation');
  case 7
    s.data.c = nf_PDNS(x);       
    fprintf('Flip-Neimark Sacker bifurcation (period = %e, parameters = %e, %e)\n',x(length(x)-3),x(length(x)-2),x(length(x)-1));
    fprintf('(p11,p22,theta,delta,sign(l1))=(%d, %d, %d, %d, %d)\n',s.data.c); 
    s.msg= sprintf('Flip-Neimark Sacker bifurcation');
    s.data.multipliers = lds.multipliers;    
  case 8
    s.data.c = nf_NSNS(x);       
    fprintf('Neimark Sacker-Neimark Sacker bifurcation (period = %e, parameters = %e, %e)\n',x(length(x)-3),x(length(x)-2),x(length(x)-1));
    fprintf('(p11,p22,theta,delta,sign(l1))=(%d, %d, %d, %d, %d)\n',s.data.c); 
    s.msg= sprintf('Neimark Sacker-Neimark Sacker bifurcation');
end
  failed = 0;
 
%--------------------------------------------------------  
function [S,L] = singmat
  S = [ 0 8 8 8 8 8 8 8
        8 0 8 8 8 8 8 8
        8 8 0 8 8 8 8 8
        8 8 8 0 8 8 8 8
        8 8 8 8 0 8 8 8
        8 8 8 8 1 0 8 8
        8 1 8 8 8 8 0 8
        8 8 8 8 8 8 8 0];
  L = [ 'R1  '; 'R2  '; 'R3  '; 'R4  ';'LPNS';'CH  ';'PDNS';'NSNS'];

%---------------------------------------------------------
function locate(varargin)
%---------------------------------------------------------
function varargout = init(varargin)

  WorkspaceInit(varargin{1:2});
  % all done succesfully
  varargout{1} = 0;
%---------------------------------------------------------
function done

  WorkspaceDone;
  
%---------------------------------------------------------
function [res,x,v] = adapt(x,v)
global lds
% calculate phi and psi for next point
[x1,p1,T,k] = rearr(x);
J = BVP_NSjac('BVP_NS_jac',x1,p1,T,k,1,1);
b = zeros((2*lds.tps-1)*lds.nphase+2,2);
b(end-1:end,:) = eye(2);  
s = J'\b;
lds.NS_new_psi = orth(s(1:end-2,:));
s = J\b;
lds.NS_new_phi = orth( s(1:end-2,:));
lds.NS_phi0 = lds.NS_new_phi(:,1);
lds.NS_phi1 = lds.NS_new_phi(:,2);
lds.NS_psi0 = lds.NS_new_psi(:,1);
lds.NS_psi1 = lds.NS_new_psi(:,2);
lds.NS1_psi = lds.NS1_new_psi(1:lds.ncoords);
lds.NS2_psi = lds.NS2_new_psi(1:lds.ncoords);
lds.NS1_phi = lds.NS1_new_phi(1:lds.ncoords);
lds.NS2_phi = lds.NS2_new_phi(1:lds.ncoords);

[x,v] = adapt_mesh(x,v);

% %compute indices
A = NS_BVP_jac('BVP_LC_jac_f','BVP_LC_jac_bc','BVP_LC_jac_ic',x1,p1,T,3,2);
[Q,R] = qr(A');
clear R
b1 = lds.NS_new_phi';
b2 = lds.NS_new_psi';
b21 = b2(:,end-lds.nphase+1:end);

% calculate g'
ups = reshape(x1,lds.nphase,lds.tps);
p   = num2cell(p1);
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
cb1 = [];
t  = lds.nphase:((lds.ncol+2)*lds.nphase-1);
kr1 = fix(t/lds.nphase);
kr2 = rem(t,lds.nphase)+1;
gx = zeros(4,lds.ncoords+3);gk=[];
for ns = 1:2
      range3 = lds.cols_p1_coords;
      b11    = cb1  ;
      range0 = lds.cols_p1;
      for tstpt = lds.tsts
          xp  = ups(:,range0)*lds.wt;
          cb1 = b1(:,range2)';
          cb2 = b2(:,range1);
          range = lds.phases;
          for c = lds.cols
              xt    = xp(:,c);
              sysj  = cjac(lds.func,lds.Jacobian,xt,p,lds.ActiveParams);
              sysh  = chess(lds.func,lds.Jacobian,lds.Hessians,xt,p,lds.ActiveParams);
              syshp = chessp(lds.func,lds.Jacobian,lds.HessiansP,xt,p,lds.ActiveParams);
              wtk   = lds.wt(kr1,c(ones(1,lds.nphase)))';
              for d = lds.phases
                  sh1(:,d) = (wtk.*sysh(:,kr2,d))*cb1(:,1);
                  sh2(:,d) = (wtk.*sysh(:,kr2,d))*cb1(:,2);
              end      
              t11 = T* wtk.*sh1(:,kr2);
              t21 = T* wtk.*sh2(:,kr2);
              t12 = (wtk.*sysj(:,kr2))*cb1(:,1);
              t22 = (wtk.*sysj(:,kr2))*cb1(:,2);
              t13 = T* wtk.*syshp(:,kr2,1)* cb1(:,1);
              t23 = T* wtk.*syshp(:,kr2,1)* cb1(:,2);
              t14 = T* wtk.*syshp(:,kr2,2)* cb1(:,1);
              t24 = T* wtk.*syshp(:,kr2,2)* cb1(:,2);
              syshess1(range,:) = [t11 t12 t13 t14];      
              syshess2(range,:) = [t21 t22 t23 t24];      
              range = range + lds.nphase;    
          end
          gx(1,[range3 lds.ncoords+(1:3)]) = gx(1,[range3 lds.ncoords+(1:3)]) + cb2(1,:)*syshess1;
          gx(2,[range3 lds.ncoords+(1:3)]) = gx(2,[range3 lds.ncoords+(1:3)]) + cb2(1,:)*syshess2;
          gx(3,[range3 lds.ncoords+(1:3)]) = gx(3,[range3 lds.ncoords+(1:3)]) + cb2(2,:)*syshess1;
          gx(4,[range3 lds.ncoords+(1:3)]) = gx(4,[range3 lds.ncoords+(1:3)]) + cb2(2,:)*syshess2;
          range0 = range0 + lds.ncol;
          range1 = range1 + lds.ncol_coord;
          range2 = range2 + lds.ncol_coord;
          range3 = range3 + lds.ncol_coord;
      end
end  
gk(1,1) = 2*b21(1,:)*b11(end-lds.nphase+1:end,1);
gk(2,1) = 2*b21(1,:)*b11(end-lds.nphase+1:end,2);
gk(3,1) = 2*b21(2,:)*b11(end-lds.nphase+1:end,1);
gk(4,1) = 2*b21(2,:)*b11(end-lds.nphase+1:end,2);
B = [A ; gx gk]*Q; clear A;
Jres = B(2+lds.ncoords:end,2+lds.ncoords:end)';
[Q,R,E] = qr(full(Jres));
index = [1 1;1 2;2 1;2 2];
[I,J] = find(E(:,1:2));
lds.index1 = index(I(J(1)),:);
lds.index2 = index(I(J(2)),:);
res = 1;


%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------

function [x,p,T,k] = rearr(x0)
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
k = x0(end);
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
function jac = BVP_jac(BVP_func,x,p,T,pars,nc)
global lds 
 
p2 = num2cell(p);
jac = feval(BVP_func,lds.func,x,p,T,pars,nc,lds,p2,lds.Jacobian,lds.ActiveParams,lds.JacobianP); 
jac = sparse(full(jac));

% -------------------------------------------------------------

function jac = BVP_NSjac(BVP_func,x,p,T,k,pars,nc)
global lds 
 
p2 = num2cell(p);
jac = feval(BVP_func,lds.func,x,p,T,k,pars,nc,lds,p2,lds.Jacobian,lds.ActiveParams,lds.JacobianP); 
jac = sparse(full(jac));


%---------------------------------------------------------------------
function [jac1,jac2] = sBVP_jac(BVP_jac_f,BVP_jac_bc,x,p,T,theta,pars,nc)
global lds 

ups = reshape(x,lds.nphase,lds.tps);
p = num2cell(p);

jac1 = spalloc(lds.ncoords,lds.ncoords,(lds.ncol+4)*lds.nphase);
%jac = zeros(lds.ncoords+1,lds.ncoords+length(p)-1);

range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;
eitheta = 1i*theta*lds.pwwt;
for jk = lds.tsts
  % value of polynomial on each collocation point
  xp = ups(:,range0)*lds.wt;

  % evaluate part of Jacobian matrix
  jac2(range1,range2) = feval(BVP_jac_f,lds.func,xp,p,T,jk);
  jac1(range1,range2) = jac2(range1,range2)+eitheta;

  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end

% boundary conditions
range = (lds.tps-1)*lds.nphase + (lds.phases);
jac1(range,[lds.phases range]) = feval(BVP_jac_bc);
jac2(range,[lds.phases range]) = feval(BVP_jac_bc);

% ---------------------------------------------------------------


function WorkspaceInit(x,v)
global cds lds

lds.cols_p1 = 1:(lds.ncol+1);
lds.cols_p1_coords = 1:(lds.ncol+1)*lds.nphase;
lds.ncol_coord = lds.ncol*lds.nphase;
lds.col_coords = 1:lds.ncol*lds.nphase;
lds.coords = 1:lds.ncoords;
lds.pars = lds.ncoords+(1:4);
lds.tsts = 1:lds.ntst;
lds.cols = 1:lds.ncol;
lds.phases = 1:lds.nphase;
lds.ntstcol = lds.ntst*lds.ncol;

lds.idxmat = reshape(fix((1:((lds.ncol+1)*lds.ntst))/(1+1/lds.ncol))+1,lds.ncol+1,lds.ntst);
lds.dt = lds.msh(lds.tsts+1)-lds.msh(lds.tsts);

lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);

[lds.NS1_phi,lds.NS2_phi,lds.NS1_psi,lds.NS2_psi] = initborders(x,v);
lds.NS1_new_phi = lds.NS1_phi;
lds.NS2_new_phi = lds.NS2_phi;
lds.NS1_new_psi = lds.NS1_psi;
lds.NS2_new_psi = lds.NS2_psi;

lds.NS_new_phi = [lds.NS_phi0 lds.NS_phi1];
lds.NS_new_psi = [lds.NS_psi0 lds.NS_psi1];

lds.NS_switch = 0;
lds.NS1_switch = 0;
lds.NS2_switch = 0;

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
function [q1,q2,p1,p2]=initborders(x,v)

global lds
[xt,p,T,k] = rearr(x);     
theta = acos(k);
[Jjac,Jjac2] = sBVP_jac('bordBVP_PD_jac_f','BVP_LC1_jac_bc',xt,p,T,theta,1,1);
%compute borders Jjac 
jace=[Jjac rand(lds.ncoords,1);rand(1,lds.ncoords) 0];
b=zeros(lds.ncoords+1,1);b(lds.ncoords+1,1)=1;
q=jace\b;q1=q(1:lds.ncoords,1);q1=q1/norm(q1);
p=jace'\b;p1=p(1:lds.ncoords,1);p1=p1/norm(p1);
jace=[Jjac2 rand(lds.ncoords,1);rand(1,lds.ncoords) 0];
q=jace\b;q2=q(1:lds.ncoords,1);q2=q2/norm(q2);
p=jace'\b;p2=p(1:lds.ncoords,1);p2=p2/norm(p2);

%SD:continues period doubling bifurcation of odefile using minimal extended system
% ---------------------------------------------------------------

