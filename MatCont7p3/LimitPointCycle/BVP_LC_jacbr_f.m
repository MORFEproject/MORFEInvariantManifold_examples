% -------------------------------------------------------------
% functions defining the jacobian of the BVP for limitcycles
% -------------------------------------------------------------

function f = BVP_LC_jacbr_f(odefile,xp,p,T,tp)
global lds
sysjacp = zeros(lds.ncol_coord,length(lds.BranchParams));
range = lds.phases;
% xp:value of polynomial on each collocation point
for c=lds.cols
  xt = xp(:,c);
  sysjacp(range,:) = cjacp(lds.func,lds.JacobianP,xt,p,lds.ActiveParams);
  range = range + lds.nphase;
end
f = -T*sysjacp;
