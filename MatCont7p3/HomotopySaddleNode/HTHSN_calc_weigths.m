function HTHSN_calc_weigths()
% calculate weights
global HTHSNds

HTHSNds.wi = nc_weight(HTHSNds.ncol)';
HTHSNds.wt = [];
HTHSNds.wpvec = [];
zm = gl_pos(HTHSNds.ncol)/2 + 0.5;
xm = (0:HTHSNds.ncol)/HTHSNds.ncol;

ncp1 = HTHSNds.ncol+1;
for j=1:ncp1
  xmi = xm(setdiff(1:ncp1,j));
  denom = prod(xm(j)-xmi);
  HTHSNds.wt(j,:) = prod(repmat(zm',HTHSNds.ncol,1)-repmat(xmi',1,HTHSNds.ncol))/denom;
  s = zeros(1,HTHSNds.ncol);
  for l=setdiff(1:ncp1,j)
    xmil = xm(setdiff(1:ncp1,[j l]));
    rep = repmat(zm',HTHSNds.ncol-1,1)-repmat(xmil',1,HTHSNds.ncol);
    if HTHSNds.ncol > 2
      s = s+prod(rep);
    else
      s = s+rep;
    end
  end
  HTHSNds.wpvec(j,:) = s/denom;
end
