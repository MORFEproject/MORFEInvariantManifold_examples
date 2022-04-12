% -------------------------------------------------------------
% functions defining the jacobian of the BVP for fold(LPC) bifurcations
% -------------------------------------------------------------

% boundary condition
function bc = BVP_LPC_bc1()
global lds
bc = [eye(lds.nphase) -eye(lds.nphase) zeros(lds.nphase,1) lds.LPC_psi(lds.ncol_coord*lds.ntst+lds.phases)'];
