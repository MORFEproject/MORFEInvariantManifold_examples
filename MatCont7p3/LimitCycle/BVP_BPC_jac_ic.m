% -------------------------------------------------------------
% functions defining the jacobian of the BVP for flip (PD) bifurcations
% -------------------------------------------------------------

% integral contraint
function ic = BVP_BPC_jac_ic()
global lds
ic=lds.BP_phi(lds.coords);
