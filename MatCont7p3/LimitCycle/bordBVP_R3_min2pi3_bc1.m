% -------------------------------------------------------------
% functions defining the jacobian of the BVP for fold(LPC) bifurcations
% -------------------------------------------------------------

% boundary condition
function bc = bordBVP_R3_min2pi3_bc1()
global lds
bc = [eye(lds.nphase) -exp(-i*2*pi/3)*eye(lds.nphase) ];
