% -------------------------------------------------------------
% functions defining the jacobian of the BVP for fold(LPC) bifurcations
% -------------------------------------------------------------

% function
function f = bordBVP_LPC_f(odefile,xp,p,T,tp)
global lds
wploc = lds.wp/lds.dt(tp);
range = lds.phases;
for c=lds.cols
  sysjac(range,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',cjac(lds.func,lds.Jacobian,xp(:,c),p,lds.ActiveParams));
  f1(range,:) = feval(odefile, 0, xp(:,c), p{:});
  range = range + lds.nphase;
end
f = [wploc-T*sysjac    -f1 ];
