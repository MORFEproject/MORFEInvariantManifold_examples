% -------------------------------------------------------------
% functions defining the jacobian of the BVP for fold (LPC) bifurcations
% -------------------------------------------------------------

% integral contraint
function bc = BVP_LPC_bc2(odefile,xp,p)
global lds
bc = zeros(1,lds.ncoords);
range0 = lds.cols_p1;
range2 = lds.cols_p1_coords;
d = zeros(lds.nphase,lds.ncol+1);
for j = lds.tsts
  xpp = xp(:,range0);
  for c = lds.cols_p1      
      d(lds.phases,c) = feval(odefile,0,xpp(:,c), p{:});
  end
  f1 = reshape((lds.dt(j)*(d.*lds.pwi))',1,lds.nphase*(lds.ncol+1));
  bc(range2) = bc(range2)+f1(lds.cols_p1_coords);
  range2 = range2 + lds.ncol_coord;
  range0 = range0 + lds.ncol;
end

bc=[bc lds.LPC_psi(end)];