% -------------------------------------------------------------
% functions defining the jacobian of the BVP for limitcycles
% -------------------------------------------------------------

% boundary condition
function bc = BVP_LC1_jac_bc()
global lds
bc = [eye(lds.nphase) -eye(lds.nphase)];
