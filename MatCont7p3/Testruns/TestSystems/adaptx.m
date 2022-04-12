function out = adaptx
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
function dydt = fun_eval(t,kmrgd,alpha,beta)
dydt=[kmrgd(2);
kmrgd(3);
-alpha*kmrgd(3)-beta*kmrgd(2)-kmrgd(1)+kmrgd(1)^2;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
y0=[0,0,0];
handles = feval(adaptx);
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,alpha,beta)
jac=[[0,1,0];[0,0,1];[-1+2*kmrgd(1),-beta,-alpha]];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,alpha,beta)
jacp=[[0,0];[0,0];[-kmrgd(3),-kmrgd(2)]];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,alpha,beta)
hess1=[[0,0,0];[0,0,0];[2,0,0]];
hess2=[[0,0,0];[0,0,0];[0,0,0]];
hess3=[[0,0,0];[0,0,0];[0,0,0]];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,alpha,beta)
hessp1=[[0,0,0];[0,0,0];[0,0,-1]];
hessp2=[[0,0,0];[0,0,0];[0,-1,0]];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,alpha,beta)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,alpha,beta)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,alpha,beta)
