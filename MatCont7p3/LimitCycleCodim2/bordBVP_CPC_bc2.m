% -------------------------------------------------------------
% functions defining the jacobian of the BVP for fold (LPC) bifurcations
% -------------------------------------------------------------

% integral contraint
function ic = bordBVP_CPC_bc2(odefile,xp,p)
global lds

range1 = lds.cols_p1;
range2 = lds.cols;
range3 = lds.cols_p1_coords;
ic = zeros(1,lds.ncoords);
wt = lds.wt';
for j = lds.tsts
  xpp = xp(:,range1)*lds.wt;
  for c = lds.cols
     xt = xpp(:,c);
     fpsC(lds.phases,c) = feval(odefile,0,xt, p{:});     
  end
  fpsW1(:,range2) = (fpsC(lds.phases,lds.cols).*repmat(gl_weight(lds.ncol)',lds.nphase,1)/2)*lds.dt(j);
  p1 = fpsW1(:,range2)*wt;
  ic(range3) = ic(range3)+reshape(p1,1,(lds.nphase*(lds.ncol+1)));
  range1 = range1 + lds.ncol;
  range2 = range2 + lds.ncol;
  range3 = range3 + lds.ncol_coord;
end