function out = bratu
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];
out{10}= @userf1;

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,a)
dydt=[-2*kmrgd(1)+kmrgd(2)+a*exp(kmrgd(1));;
kmrgd(1)-2*kmrgd(2)+a*exp(kmrgd(2));;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
handles = feval(bratu);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,a)
jac=[ a*exp(kmrgd(1)) - 2 , 1 ; 1 , a*exp(kmrgd(2)) - 2 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,a)
jacp=[ exp(kmrgd(1)) ; exp(kmrgd(2)) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,a)
hess1=[ a*exp(kmrgd(1)) , 0 ; 0 , 0 ];
hess2=[ 0 , 0 ; 0 , a*exp(kmrgd(2)) ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,a)
hessp1=[ exp(kmrgd(1)) , 0 ; 0 , exp(kmrgd(2)) ];
hessp(:,:,1) =hessp1;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,a)
tens31=[ a*exp(kmrgd(1)) , 0 ; 0 , 0 ];
tens32=[ 0 , 0 ; 0 , 0 ];
tens33=[ 0 , 0 ; 0 , 0 ];
tens34=[ 0 , 0 ; 0 , a*exp(kmrgd(2)) ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,2,1) =tens33;
tens3(:,:,2,2) =tens34;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,a)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,a)
function userfun1=userf1(t,kmrgd,a)
	userfun1=a-0.2;
