% -------------------------------------------------------------
% functions defining the jacobian of the BVP for fold(LPC) bifurcations
% -------------------------------------------------------------

% boundary condition
function bc = bordBVP_NS_bc1(k)
global lds
bc = [eye(lds.nphase) -2*k*eye(lds.nphase) eye(lds.nphase)];
