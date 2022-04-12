% -------------------------------------------------------------
% functions defining the jacobian of the BVP for fold(LPC) bifurcations
% -------------------------------------------------------------

% boundary condition
function bc = bordBVP_LPC_bc1()
global lds
bc = [eye(lds.nphase) -eye(lds.nphase) ];
