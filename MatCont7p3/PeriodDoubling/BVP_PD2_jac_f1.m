% -------------------------------------------------------------
% functions defining the jacobian of the BVP for limitcycles
% -------------------------------------------------------------

% function
function f = BVP_LC_jac_f1(odefile,xp,yp,p,T,tp)
global lds
sysjacp = zeros(lds.ncol_coord,length(lds.ActiveParams));
wploc = lds.wp/lds.dt(tp);
range = lds.phases;
for c=lds.cols
  xt = xp(:,c);
  yt = yp(:,c);
  frhs(range) = cjac(lds.func,lds.Jacobian,xt,p,lds.ActiveParams)*yt;
  hess = chess(lds.func,lds.Jacobian,lds.Hessians,xt,p,lds.ActiveParams);
  for i=1:size(hess,3)
    h(:,i) = hess(:,:,i)*yt;
  end
  sysjac(range,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',h);
  jacp = chessp(lds.func,lds.Jacobian,lds.HessiansP,xt,p,lds.ActiveParams);
  for i=1:size(jacp,3)
    sysjp(:,i) = jacp(:,:,i)*yt;
  end
  sysjacp(range,:) = sysjp(:,:);
  range = range + lds.nphase;
end
f = [-T*sysjac    -frhs'    -T*sysjacp];
