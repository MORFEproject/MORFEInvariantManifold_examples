% -------------------------------------------------------------
% functions defining the BVP for limitcycles
% -------------------------------------------------------------

function ic = BVP_PD2_ic(ups)
global lds
ficd = dot(ups,lds.vpoldp);
ficdmat = ficd(lds.idxmat);
ic = sum(lds.dt.*(lds.wi*ficdmat))-1;
