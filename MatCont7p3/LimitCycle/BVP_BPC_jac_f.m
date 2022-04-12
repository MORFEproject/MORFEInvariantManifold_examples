% -------------------------------------------------------------
% functions defining the jacobian of the BVP for BPC bifurcations
% -------------------------------------------------------------

% function
function f = BVP_BPC_jac_f(odefile,xp,p,T,tp)
global lds
wploc = lds.wp/lds.dt(tp);
range = lds.phases;
for c=lds.cols
  sysjac(range,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',cjac(lds.func,lds.Jacobian,xp(:,c),p,lds.ActiveParams));
  range = range + lds.nphase;
end
f = [wploc-T*sysjac    lds.BP_psi((tp-1)*lds.ncol_coord + lds.col_coords)'];
