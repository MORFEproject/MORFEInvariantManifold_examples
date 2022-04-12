function out = cstr
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,lambda,beta,gamma,alpha4)
alpha1=10-9*beta+gamma;
alpha2=10-9*beta;
alpha3=-0.9+0.4*beta;
dydt=[alpha3-(1+lambda)*kmrgd(1)+lambda*alpha1/(1+lambda*alpha2*exp(-alpha4*kmrgd(1)/(1+kmrgd(1))));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(cstr);
y0=[0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,lambda,beta,gamma,alpha4)
jac=[ - lambda - (lambda^2*(9*beta - 10)*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) - 1 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,lambda,beta,gamma,alpha4)
jacp=[ (lambda*(9*beta - 10)*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) - (gamma - 9*beta + 10)/((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1) - kmrgd(1) , (9*lambda)/((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1) + (9*lambda^2*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) + 2/5 , -lambda/((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1) , -(lambda^2*kmrgd(1)*(9*beta - 10)*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2*(kmrgd(1) + 1)) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,lambda,beta,gamma,alpha4)
hess1=[ (lambda^2*(9*beta - 10)*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)^2*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) - (2*lambda^3*(9*beta - 10)^2*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)^2*(gamma - 9*beta + 10))/(exp((2*alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^3) + (lambda^2*(9*beta - 10)*((2*alpha4)/(kmrgd(1) + 1)^2 - (2*alpha4*kmrgd(1))/(kmrgd(1) + 1)^3)*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) ];
hess(:,:,1) =hess1;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,lambda,beta,gamma,alpha4)
hessp1=[ (2*lambda^2*(9*beta - 10)^2*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)*(gamma - 9*beta + 10))/(exp((2*alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^3) - (2*lambda*(9*beta - 10)*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) - 1 ];
hessp2=[ (9*lambda^2*(9*beta - 10)*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) - (9*lambda^2*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) + (18*lambda^3*(9*beta - 10)*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)*(gamma - 9*beta + 10))/(exp((2*alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^3) ];
hessp3=[ -(lambda^2*(9*beta - 10)*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) ];
hessp4=[ (lambda^2*(kmrgd(1)/(kmrgd(1) + 1)^2 - 1/(kmrgd(1) + 1))*(9*beta - 10)*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) - (2*lambda^3*kmrgd(1)*(9*beta - 10)^2*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)*(gamma - 9*beta + 10))/(exp((2*alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^3*(kmrgd(1) + 1)) + (lambda^2*kmrgd(1)*(9*beta - 10)*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2*(kmrgd(1) + 1)) ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,lambda,beta,gamma,alpha4)
tens31=[ (2*lambda^3*(9*beta - 10)^2*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)^3*(gamma - 9*beta + 10))/(exp((2*alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^3) - (lambda^2*(9*beta - 10)*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)^3*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) - (6*lambda^4*(9*beta - 10)^3*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)^3*(gamma - 9*beta + 10))/(exp((3*alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^4) - (lambda^2*(9*beta - 10)*((6*alpha4)/(kmrgd(1) + 1)^3 - (6*alpha4*kmrgd(1))/(kmrgd(1) + 1)^4)*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) + (6*lambda^3*(9*beta - 10)^2*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)*((2*alpha4)/(kmrgd(1) + 1)^2 - (2*alpha4*kmrgd(1))/(kmrgd(1) + 1)^3)*(gamma - 9*beta + 10))/(exp((2*alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^3) + (2*lambda^3*(9*beta - 10)^2*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)^2*((2*alpha4)/(kmrgd(1) + 1) - (2*alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)*(gamma - 9*beta + 10))/(exp((2*alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^3) - (3*lambda^2*(9*beta - 10)*(alpha4/(kmrgd(1) + 1) - (alpha4*kmrgd(1))/(kmrgd(1) + 1)^2)*((2*alpha4)/(kmrgd(1) + 1)^2 - (2*alpha4*kmrgd(1))/(kmrgd(1) + 1)^3)*(gamma - 9*beta + 10))/(exp((alpha4*kmrgd(1))/(kmrgd(1) + 1))*((lambda*(9*beta - 10))/exp((alpha4*kmrgd(1))/(kmrgd(1) + 1)) - 1)^2) ];
tens3(:,:,1,1) =tens31;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,lambda,beta,gamma,alpha4)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,lambda,beta,gamma,alpha4)
