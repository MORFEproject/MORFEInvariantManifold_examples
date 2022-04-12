% -------------------------------------------------------------
% functions defining the jacobian of the BVP for fold (LPC) bifurcations
% -------------------------------------------------------------

% integral contraint

function J = BVP_jac_LPC_ic( ups, p, cv)

global lds
t   = lds.nphase:((lds.ncol+1)*lds.nphase-1);
kr  = fix(t/lds.nphase);
wi  = lds.wi(1,kr) ;
wit = lds.wi;
wit = wit(1,kr);
t   = lds.nphase*lds.ncol:((lds.ncol)*lds.nphase*(lds.ntst)-1);
kr  = rem(t,lds.nphase*lds.ncol) + 1;
tend = lds.wi(:,end);
reswi = cat(2,wi,wit(1,kr),tend(:,ones(1,lds.nphase)));
t   = lds.ncol*lds.nphase:((lds.ncol*lds.nphase*(lds.ntst+1))-1);
kr  = fix(t/(lds.ncol*lds.nphase));
ind = [kr,kr(end)*ones(1,lds.nphase)];
res1 = reswi.*lds.dt(1,ind);
res2 = zeros(1,lds.tps);
res2(1+lds.ncol:lds.ncol:end-1) = lds.wi(end)*lds.dt(1:lds.ntst-1);
t   = lds.nphase:((lds.tps+1)*(lds.nphase)-1);
kr  = fix(t/(lds.nphase));
res = res1 + res2(1,kr);
range = lds.phases;
for d = 1:lds.tps
      jac  = cjac(lds.func,lds.Jacobian,ups(:,d),p,lds.ActiveParams);
      jacp = cjacp(lds.func,lds.JacobianP,ups(:,d),p,lds.ActiveParams);
      jacp1ic(:,range) = jacp(:,1)';
      jacp2ic(:,range) = jacp(:,2)';
      jacic(:,range)   = (jac'*cv(range,:))';
      range = range + lds.phases;
end
jacic   = res .* jacic;
jacp1ic = (res .* jacp1ic)*cv;
jacp2ic = (res .* jacp2ic)*cv;
J=[jacic jacp1ic jacp2ic];
