function out = branchpoint
%
% Branch point curve definition file for a problem in odefile
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
global  bpds
  [x,p] = rearr(arg); p = num2cell(p);
  Bord=[cjac(bpds.func,bpds.Jacobian,x,p,bpds.ActiveParams) cjacbr(bpds.func,bpds.JacobianP,x,p,bpds.ActiveParams,bpds.BranchParam) bpds.borders.w;bpds.borders.v' zeros(2,1)];
  bunit=[zeros(bpds.nphase,2);eye(2)];
  vext=Bord\bunit;
  func = [feval(bpds.func, 0, x, p{:}); vext(end,1);vext(end,2)];

%------------------------------------------------------
function jac = jacobian(varargin)
global  bpds
  nap = length(bpds.ActiveParams);
  xo = varargin{1}; [x,p] = rearr(xo); p = num2cell(p);
  jac = cjac(bpds.func,bpds.Jacobian,x,p,bpds.ActiveParams);
  Bord = [jac cjacbr(bpds.func,bpds.JacobianP,x,p,bpds.ActiveParams,bpds.BranchParam) bpds.borders.w;bpds.borders.v' zeros(2,1)];
  bunit = [zeros(bpds.nphase,2);eye(2)];
  vext = Bord\bunit;
  bunit=[zeros(bpds.nphase+1,1);1];
  wext = Bord'\bunit;
  jac = [jac cjacp(bpds.func,bpds.JacobianP,x,p,bpds.ActiveParams)];
  hess  = chess(bpds.func,bpds.Jacobian,bpds.Hessians,x,p,bpds.ActiveParams);
  hessp = chessp(bpds.func,bpds.Jacobian,bpds.HessiansP,x,p,bpds.ActiveParams);
  hessbr= chessbr(bpds.func,bpds.Jacobian,bpds.HessiansP,x,p,bpds.ActiveParams,bpds.BranchParam);
  for i=1:bpds.nphase
     jac(bpds.nphase+1,i)=-wext(1:bpds.nphase)'*hess(:,:,i)*vext(1:bpds.nphase,1)-wext(1:bpds.nphase)'*hessbr(:,i)*vext(end-1,1);
     jac(bpds.nphase+2,i)=-wext(1:bpds.nphase)'*hess(:,:,i)*vext(1:bpds.nphase,2)-wext(1:bpds.nphase)'*hessbr(:,i)*vext(end-1,2);
  end
  hesspp = chesspbr(bpds.func,bpds.JacobianP,bpds.HessiansP,x,p,bpds.ActiveParams,bpds.BranchParam);
  for i=1:nap
      jac(bpds.nphase+1,bpds.nphase+i)=-wext(1:bpds.nphase)'*hessp(:,:,i)*vext(1:bpds.nphase,1)-wext(1:bpds.nphase)'*hesspp(:,:,i)*vext(end-1,1);
      jac(bpds.nphase+2,bpds.nphase+i)=-wext(1:bpds.nphase)'*hessp(:,:,i)*vext(1:bpds.nphase,2)-wext(1:bpds.nphase)'*hesspp(:,:,i)*vext(end-1,2);     
  end
  
%-------------------------------------------------------
function hessians()
  
%-------------------------------------------------------  
function varargout = defaultprocessor(varargin)
global cds bpds
  if nargin > 2
    s = varargin{3};
    varargout{3} = s;
  end

 % compute eigenvalues?
  if (cds.options.Eigenvalues==1)
      xo = varargin{1}; [x,p] = rearr(xo); p = num2cell(p);
      d = eig(cjac(bpds.func,bpds.Jacobian,x,p,bpds.ActiveParams));
      [Y,I] = sort(real(d));
      varargout{2} = d(I);
  else
      varargout{2} = nan;
  end  

  % all done succesfully
  varargout{1} = 0;

%----------------------------------------------------------
function option = options
global bpds
  option = contset;

  % Check for symbolic derivatives in odefile
  symord = 0; 
  if ~isempty(bpds.Jacobian), symord = 1; end
  if ~isempty(bpds.Hessians), symord = 2; end
  if ~isempty(bpds.Der3), symord = 3; end
  if ~isempty(bpds.Der4), symord = 4; end
  if ~isempty(bpds.Der5), symord = 5; end
  option = contset(option, 'SymDerivative', symord);

  symordp = 0;
  if ~isempty(bpds.HessiansP)
      symordp = 2;
  elseif ~isempty(bpds.JacobianP)  
      symordp = 1; 
  end
  option = contset(option,'SymDerivativeP',symordp);

  option = contset(option, 'Workspace', 1);
  option = contset(option, 'Locators', [0 0 0]);
  option = contset(option, 'Singularities', 0);
  
%-----------------------------------------------------
function testf(varargin)

%-----------------------------------------------------
function [out, failed] = userf(userinf, id, x, v)
global  bpds
dim =size(id,2);
failed = [];
out = zeros(1,dim);
for i=1:dim
  lastwarn('');
  [x0,p] = rearr(x); p = num2cell(p);
  if (userinf(i).state==1)
      out(i)=feval(bpds.user{id(i)},0,x0,p{:});
  end
  if ~isempty(lastwarn)
    failed = [failed i];
  end
end
%---------------------------------------------------------
function process(varargin)
%---------------------------------------------------------
function singmat(varargin)
%---------------------------------------------------------
function locate(varargin)
%---------------------------------------------------------
function varargout = init(varargin)
  x = varargin{1};
  v = varargin{2};
  WorkspaceInit(x,v);

  % all done succesfully
  varargout{1} = 0;
%---------------------------------------------------------
function done
  WorkspaceDone;

% -------------------------------------------------------
function [res,x,v] = adapt(x,v)
global bpds
[x1,p] =rearr(x); p = num2cell(p);
Bord=[cjac(bpds.func,bpds.Jacobian,x1,p,bpds.ActiveParams) cjacbr(bpds.func,bpds.JacobianP,x1,p,bpds.ActiveParams,bpds.BranchParam) bpds.borders.w;bpds.borders.v' zeros(2,1)];
bunit=[zeros(bpds.nphase,2);eye(2)];
vext=Bord\bunit;
bunit=[zeros(bpds.nphase+1,1);1];
wext = Bord'\bunit;
bpds.borders.w=wext(1:bpds.nphase);
bpds.borders.v=vext(1:bpds.nphase+1,:)/norm(vext(1:bpds.nphase+1,:));
res = []; % no re-evaluations needed

  
%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------

function [x,p] = rearr(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x) and parameters (p)
global bpds
p = bpds.P0;
p(bpds.ActiveParams) = x0((bpds.nphase+1):end);
x = x0(1:bpds.nphase);

% ---------------------------------------------------------

function WorkspaceInit(x,v)

% ------------------------------------------------------

function WorkspaceDone

%SD:continues branchpoint of odefile
