function out = bratu2(t,coordinates,flag,a,b,c)
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,a,b,c)
dydt=[kmrgd(2)-2*kmrgd(1)+a*exp(kmrgd(1))+b+c;
kmrgd(1)-2*kmrgd(2)+a*exp(kmrgd(2))+b;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(bratu2);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,a,b,c)
jac=[[-2+a*exp(kmrgd(1)),1];[1,-2+a*exp(kmrgd(2))]];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,a,b,c)
jacp=[[exp(kmrgd(1)),1,1];[exp(kmrgd(2)),1,0]];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,a,b,c)
hess1=[[a*exp(kmrgd(1)),0];[0,0]];
hess2=[[0,0];[0,a*exp(kmrgd(2))]];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,a,b,c)
hessp1=[[exp(kmrgd(1)),0];[0,exp(kmrgd(2))]];
hessp2=[[0,0];[0,0]];
hessp3=[[0,0];[0,0]];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,a,b,c)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,a,b,c)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,a,b,c)
