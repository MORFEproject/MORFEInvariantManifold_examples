function Hom_calc_weights()
% calculate weights

global homds

homds.wi = nc_weight(homds.ncol)';
homds.wt = [];
homds.wpvec = [];
zm = gl_pos(homds.ncol)/2 + 0.5;
xm = (0:homds.ncol)/homds.ncol;

ncp1 = homds.ncol+1;
for j=1:ncp1
  xmi = xm(setdiff(1:ncp1,j));
  denom = prod(xm(j)-xmi);
  homds.wt(j,:) = prod(repmat(zm',homds.ncol,1)-repmat(xmi',1,homds.ncol))/denom;
  s = zeros(1,homds.ncol);
  for l=setdiff(1:ncp1,j)
    xmil = xm(setdiff(1:ncp1,[j l]));
    rep = repmat(zm',homds.ncol-1,1)-repmat(xmil',1,homds.ncol);
    if homds.ncol > 2
      s = s+prod(rep);
    else
      s = s+rep;
    end
  end
  homds.wpvec(j,:) = s/denom;
end
