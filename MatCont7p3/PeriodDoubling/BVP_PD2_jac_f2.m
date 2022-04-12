% -------------------------------------------------------------
% functions defining the jacobian of the BVP for limitcycles
% -------------------------------------------------------------

% function
function f = BVP_PD2_jac_f2(odefile,xp,yp,p,T,tp)
global lds
wploc = lds.wp/lds.dt(tp);
range = lds.phases;
for c=lds.cols
  xt = xp(:,c);
  yt = yp(:,c);
  sysjac(range,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',cjac(lds.func,lds.Jacobian,xt,p,lds.ActiveParams));
  range = range + lds.nphase;
end
f = [wploc-T*sysjac];
