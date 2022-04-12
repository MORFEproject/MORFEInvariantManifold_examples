% -------------------------------------------------------------
% functions defining the jacobian of the BVP for flip (PD) bifurcations
% -------------------------------------------------------------

% boundary condition
function bc = bordBVP_PD_jac_bc()
global lds
bc = [eye(lds.nphase) eye(lds.nphase) ];
