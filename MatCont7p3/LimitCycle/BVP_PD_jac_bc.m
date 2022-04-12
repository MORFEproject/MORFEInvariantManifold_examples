% -------------------------------------------------------------
% functions defining the jacobian of the BVP for flip (PD) bifurcations
% -------------------------------------------------------------

% boundary condition
function bc = BVP_PD_jac_bc()
global lds
bc = [eye(lds.nphase) eye(lds.nphase) lds.PD_psi(lds.ncol_coord*lds.ntst+lds.phases)'];
