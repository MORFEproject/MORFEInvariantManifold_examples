function HTHet_calc_weigths()
% calculate weights
global HTHetds

HTHetds.wi = nc_weight(HTHetds.ncol)';
HTHetds.wt = [];
HTHetds.wpvec = [];
zm = gl_pos(HTHetds.ncol)/2 + 0.5;
xm = (0:HTHetds.ncol)/HTHetds.ncol;

ncp1 = HTHetds.ncol+1;
for j=1:ncp1
  xmi = xm(setdiff(1:ncp1,j));
  denom = prod(xm(j)-xmi);
  HTHetds.wt(j,:) = prod(repmat(zm',HTHetds.ncol,1)-repmat(xmi',1,HTHetds.ncol))/denom;
  s = zeros(1,HTHetds.ncol);
  for l=setdiff(1:ncp1,j)
    xmil = xm(setdiff(1:ncp1,[j l]));
    rep = repmat(zm',HTHetds.ncol-1,1)-repmat(xmil',1,HTHetds.ncol);
    if HTHetds.ncol > 2
      s = s+prod(rep);
    else
      s = s+rep;
    end
  end
  HTHetds.wpvec(j,:) = s/denom;
end
