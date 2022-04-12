% -------------------------------------------------------------
% functions defining the BVP for limitcycles
% -------------------------------------------------------------

function ic = BVP_LC_ic(ups)
global lds
ficd = dot(ups,lds.upoldp);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat));
