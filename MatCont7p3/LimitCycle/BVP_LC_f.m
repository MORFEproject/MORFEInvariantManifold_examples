% -------------------------------------------------------------
% functions defining the BVP for limitcycles
% -------------------------------------------------------------

function f = BVP_LC_f(odefile,t,xp,p,T)
f = t-T*feval(odefile, 0, xp, p{:});
