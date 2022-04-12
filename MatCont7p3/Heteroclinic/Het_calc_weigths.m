function Het_calc_weigths()
% calculate weights
global hetds

hetds.wi = nc_weight(hetds.ncol)';
hetds.wt = [];
hetds.wpvec = [];
zm = gl_pos(hetds.ncol)/2 + 0.5;
xm = (0:hetds.ncol)/hetds.ncol;

ncp1 = hetds.ncol+1;
for j=1:ncp1
  xmi = xm(setdiff(1:ncp1,j));
  denom = prod(xm(j)-xmi);
  hetds.wt(j,:) = prod(repmat(zm',hetds.ncol,1)-repmat(xmi',1,hetds.ncol))/denom;
  s = zeros(1,hetds.ncol);
  for l=setdiff(1:ncp1,j)
    xmil = xm(setdiff(1:ncp1,[j l]));
    rep = repmat(zm',hetds.ncol-1,1)-repmat(xmil',1,hetds.ncol);
    if hetds.ncol > 2
      s = s+prod(rep);
    else
      s = s+rep;
    end
  end
  hetds.wpvec(j,:) = s/denom;
end
