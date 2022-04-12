% -------------------------------------------------------------
% functions defining the jacobian of the BVP for limitcycles
% -------------------------------------------------------------
function f = BVP_LC_jac_f(odefile,xp,p,T,tp)
global lds
sysjacp = zeros(lds.ncol_coord,length(lds.ActiveParams));
wploc = lds.wp/lds.dt(tp);
range = lds.phases;
% xp:value of polynomial on each collocation point
for c=lds.cols
  xt = xp(:,c);
  frhs(range) = feval(odefile, 0, xt, p{:});
  sysjac(range,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',cjac(lds.func,lds.Jacobian,xt,p,lds.ActiveParams));
  sysjacp(range,:) = cjacp(lds.func,lds.JacobianP,xt,p,lds.ActiveParams);
  range = range + lds.nphase;
end
f = [wploc-T*sysjac    -frhs'    -T*sysjacp];
