function out = MyML
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
function dydt = fun_eval(t,kmrgd,Inp,V3)
Minf=(1+tanh((kmrgd(1)+1.2)/18))/2;
Ninf=(1+tanh((kmrgd(1)-V3)/17.4))/2;
tau=1/15*cosh((kmrgd(1)-V3)/34.8);
dydt=[1/5*(Inp-2*(kmrgd(1)+60)-4*Minf*(kmrgd(1)-120)-8*kmrgd(2)*(kmrgd(1)+80));
tau*(Ninf-kmrgd(2));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(MyML);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,Inp,V3)
jac=[[-4/5-1/5*(1/9-1/9*tanh(1/18*kmrgd(1)+1/15)^2)*(kmrgd(1)-120)-2/5*tanh(1/18*kmrgd(1)+1/15)-8/5*kmrgd(2),-8/5*kmrgd(1)-128];[1/522*sinh(5/174*kmrgd(1)-5/174*V3)*(1/2+1/2*tanh(5/87*kmrgd(1)-5/87*V3)-kmrgd(2))+1/15*cosh(5/174*kmrgd(1)-5/174*V3)*(5/174-5/174*tanh(5/87*kmrgd(1)-5/87*V3)^2),-1/15*cosh(5/174*kmrgd(1)-5/174*V3)]];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,Inp,V3)
jacp=[[1/5,0];[0,-1/522*sinh(5/174*kmrgd(1)-5/174*V3)*(1/2+1/2*tanh(5/87*kmrgd(1)-5/87*V3)-kmrgd(2))+1/15*cosh(5/174*kmrgd(1)-5/174*V3)*(-5/174+5/174*tanh(5/87*kmrgd(1)-5/87*V3)^2)]];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,Inp,V3)
hess1=[[2/45*tanh(1/18*kmrgd(1)+1/15)*(1/18-1/18*tanh(1/18*kmrgd(1)+1/15)^2)*(kmrgd(1)-120)-2/45+2/45*tanh(1/18*kmrgd(1)+1/15)^2,-8/5];[5/90828*cosh(5/174*kmrgd(1)-5/174*V3)*(1/2+1/2*tanh(5/87*kmrgd(1)-5/87*V3)-kmrgd(2))+1/261*sinh(5/174*kmrgd(1)-5/174*V3)*(5/174-5/174*tanh(5/87*kmrgd(1)-5/87*V3)^2)-1/261*cosh(5/174*kmrgd(1)-5/174*V3)*tanh(5/87*kmrgd(1)-5/87*V3)*(5/87-5/87*tanh(5/87*kmrgd(1)-5/87*V3)^2),-1/522*sinh(5/174*kmrgd(1)-5/174*V3)]];
hess2=[[-8/5,0];[-1/522*sinh(5/174*kmrgd(1)-5/174*V3),0]];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,Inp,V3)
hessp1=[[0,0];[0,0]];
hessp2=[[0,0];[-5/90828*cosh(5/174*kmrgd(1)-5/174*V3)*(1/2+1/2*tanh(5/87*kmrgd(1)-5/87*V3)-kmrgd(2))+1/522*sinh(5/174*kmrgd(1)-5/174*V3)*(-5/174+5/174*tanh(5/87*kmrgd(1)-5/87*V3)^2)-1/522*sinh(5/174*kmrgd(1)-5/174*V3)*(5/174-5/174*tanh(5/87*kmrgd(1)-5/87*V3)^2)-1/261*cosh(5/174*kmrgd(1)-5/174*V3)*tanh(5/87*kmrgd(1)-5/87*V3)*(-5/87+5/87*tanh(5/87*kmrgd(1)-5/87*V3)^2),1/522*sinh(5/174*kmrgd(1)-5/174*V3)]];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,Inp,V3)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,Inp,V3)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,Inp,V3)
