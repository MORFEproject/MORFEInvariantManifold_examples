% -------------------------------------------------------------
% functions defining the jacobian of the BVP for fold(LPC) bifurcations
% -------------------------------------------------------------

% boundary condition
function bc = bordBVP_R4_minpi2_bc1()
global lds
bc = [eye(lds.nphase) -exp(-i*pi/2)*eye(lds.nphase) ];
