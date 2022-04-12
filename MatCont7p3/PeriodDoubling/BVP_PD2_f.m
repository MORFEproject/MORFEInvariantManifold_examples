% -------------------------------------------------------------
% functions defining the BVP for limitcycles
% -------------------------------------------------------------

function f = BVP_PD2_f(odefile,t,xp,yp,p,T)
global lds
f = t-T*cjac(lds.func,lds.Jacobian,xp,p,lds.ActiveParams)*yp;

