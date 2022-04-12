function out = torBPC
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];
out{10}= @epsilon0;

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,beta,gamma,r,a3,b3,nu,epsilon)
dydt=[(-(beta+nu)*kmrgd(1)+beta*kmrgd(2)-a3*kmrgd(1)^3+b3*(kmrgd(2)-kmrgd(1))^3)/r;
beta*kmrgd(1)-(beta+gamma)*kmrgd(2)-kmrgd(3)-b3*(kmrgd(2)-kmrgd(1))^3;
kmrgd(2)+epsilon;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
y0=[0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
handles = feval(torBPC);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,beta,gamma,r,a3,b3,nu,epsilon)
jac=[ -(beta + nu + 3*a3*kmrgd(1)^2 + 3*b3*(kmrgd(1) - kmrgd(2))^2)/r , (beta + 3*b3*(kmrgd(1) - kmrgd(2))^2)/r , 0 ; beta + 3*b3*(kmrgd(1) - kmrgd(2))^2 , - beta - gamma - 3*b3*(kmrgd(1) - kmrgd(2))^2 , -1 ; 0 , 1 , 0 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,beta,gamma,r,a3,b3,nu,epsilon)
jacp=[ -(kmrgd(1) - kmrgd(2))/r , 0 , (kmrgd(1)*(beta + nu) - beta*kmrgd(2) + a3*kmrgd(1)^3 + b3*(kmrgd(1) - kmrgd(2))^3)/r^2 , -kmrgd(1)^3/r , -(kmrgd(1) - kmrgd(2))^3/r , -kmrgd(1)/r , 0 ; kmrgd(1) - kmrgd(2) , -kmrgd(2) , 0 , 0 , (kmrgd(1) - kmrgd(2))^3 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 , 1 ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,beta,gamma,r,a3,b3,nu,epsilon)
hess1=[ -(6*a3*kmrgd(1) + 3*b3*(2*kmrgd(1) - 2*kmrgd(2)))/r , (3*b3*(2*kmrgd(1) - 2*kmrgd(2)))/r , 0 ; 3*b3*(2*kmrgd(1) - 2*kmrgd(2)) , -3*b3*(2*kmrgd(1) - 2*kmrgd(2)) , 0 ; 0 , 0 , 0 ];
hess2=[ (3*b3*(2*kmrgd(1) - 2*kmrgd(2)))/r , -(3*b3*(2*kmrgd(1) - 2*kmrgd(2)))/r , 0 ; -3*b3*(2*kmrgd(1) - 2*kmrgd(2)) , 3*b3*(2*kmrgd(1) - 2*kmrgd(2)) , 0 ; 0 , 0 , 0 ];
hess3=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,beta,gamma,r,a3,b3,nu,epsilon)
hessp1=[ -1/r , 1/r , 0 ; 1 , -1 , 0 ; 0 , 0 , 0 ];
hessp2=[ 0 , 0 , 0 ; 0 , -1 , 0 ; 0 , 0 , 0 ];
hessp3=[ (beta + nu + 3*a3*kmrgd(1)^2 + 3*b3*(kmrgd(1) - kmrgd(2))^2)/r^2 , -(beta + 3*b3*(kmrgd(1) - kmrgd(2))^2)/r^2 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp4=[ -(3*kmrgd(1)^2)/r , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp5=[ -(3*(kmrgd(1) - kmrgd(2))^2)/r , (3*(kmrgd(1) - kmrgd(2))^2)/r , 0 ; 3*(kmrgd(1) - kmrgd(2))^2 , -3*(kmrgd(1) - kmrgd(2))^2 , 0 ; 0 , 0 , 0 ];
hessp6=[ -1/r , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp7=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
hessp(:,:,6) =hessp6;
hessp(:,:,7) =hessp7;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,beta,gamma,r,a3,b3,nu,epsilon)
tens31=[ -(6*a3 + 6*b3)/r , (6*b3)/r , 0 ; 6*b3 , -6*b3 , 0 ; 0 , 0 , 0 ];
tens32=[ (6*b3)/r , -(6*b3)/r , 0 ; -6*b3 , 6*b3 , 0 ; 0 , 0 , 0 ];
tens33=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens34=[ (6*b3)/r , -(6*b3)/r , 0 ; -6*b3 , 6*b3 , 0 ; 0 , 0 , 0 ];
tens35=[ -(6*b3)/r , (6*b3)/r , 0 ; 6*b3 , -6*b3 , 0 ; 0 , 0 , 0 ];
tens36=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens37=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens38=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens39=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,1,3) =tens33;
tens3(:,:,2,1) =tens34;
tens3(:,:,2,2) =tens35;
tens3(:,:,2,3) =tens36;
tens3(:,:,3,1) =tens37;
tens3(:,:,3,2) =tens38;
tens3(:,:,3,3) =tens39;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,beta,gamma,r,a3,b3,nu,epsilon)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,beta,gamma,r,a3,b3,nu,epsilon)
function userfun1=epsilon0(t,kmrgd,beta,gamma,r,a3,b3,nu,epsilon)
	userfun1=epsilon;
