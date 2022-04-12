function calc_weigths()
% calculate weights

global lds

lds.wi = nc_weight(lds.ncol)';
lds.wt = [];
lds.wpvec = [];
zm = gl_pos(lds.ncol)/2 + 0.5;
xm = (0:lds.ncol)/lds.ncol;

ncp1 = lds.ncol+1;
for j=1:ncp1
  xmi = xm(setdiff(1:ncp1,j));
  denom = prod(xm(j)-xmi);
  lds.wt(j,:) = prod(repmat(zm',lds.ncol,1)-repmat(xmi',1,lds.ncol))/denom;
  s = zeros(1,lds.ncol);
  for l=setdiff(1:ncp1,j)
    xmil = xm(setdiff(1:ncp1,[j l]));
    rep = repmat(zm',lds.ncol-1,1)-repmat(xmil',1,lds.ncol);
    if lds.ncol > 2
      s = s+prod(rep);
    else
      s = s+rep;
    end
  end
  lds.wpvec(j,:) = s/denom;
end
