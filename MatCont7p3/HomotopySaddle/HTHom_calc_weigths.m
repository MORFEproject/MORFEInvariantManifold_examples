function HTHom_calc_weigths()
% calculate weights
global HTHomds

HTHomds.wi = nc_weight(HTHomds.ncol)';
HTHomds.wt = [];
HTHomds.wpvec = [];
zm = gl_pos(HTHomds.ncol)/2 + 0.5;
xm = (0:HTHomds.ncol)/HTHomds.ncol;

ncp1 = HTHomds.ncol+1;
for j=1:ncp1
  xmi = xm(setdiff(1:ncp1,j));
  denom = prod(xm(j)-xmi);
  HTHomds.wt(j,:) = prod(repmat(zm',HTHomds.ncol,1)-repmat(xmi',1,HTHomds.ncol))/denom;
  s = zeros(1,HTHomds.ncol);
  for l=setdiff(1:ncp1,j)
    xmil = xm(setdiff(1:ncp1,[j l]));
    rep = repmat(zm',HTHomds.ncol-1,1)-repmat(xmil',1,HTHomds.ncol);
    if HTHomds.ncol > 2
      s = s+prod(rep);
    else
      s = s+rep;
    end
  end
  HTHomds.wpvec(j,:) = s/denom;
end
