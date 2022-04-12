% -------------------------------------------------------------
% functions defining the jacobian of the BVP for fold(LPC) bifurcations
% -------------------------------------------------------------

% boundary condition
function bc = BVP_BPC_bc1()
global lds
bc = [eye(lds.nphase) -eye(lds.nphase) zeros(lds.nphase,1) zeros(lds.nphase,1) lds.BPC_psi(lds.ncol_coord*lds.ntst+lds.phases)'];
