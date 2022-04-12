function jac = NS_BVP_jac(BVP_jac_f,BVP_jac_bc,BVP_jac_ic,x,p,T,npar,nc)
global lds

ups = reshape(x,lds.nphase,lds.tps);
p = num2cell(p);
pars = lds.ncoords+(1:npar);
jac = spalloc(lds.ncoords+1,lds.ncoords+4,(lds.ncol+4)*lds.nphase);

range0 = lds.cols_p1;
range1 = lds.col_coords;
range2 = lds.cols_p1_coords;

for j=lds.tsts
  xp = ups(:,range0)*lds.wt;

  jac(range1,[range2 pars]) = feval(BVP_jac_f,lds.func,xp,p,T,j);

  range0 = range0 + lds.ncol;
  range1 = range1 + lds.ncol_coord;
  range2 = range2 + lds.ncol_coord;
end

% integral constraint
jac(lds.ncoords+1,[lds.coords]) = feval(BVP_jac_ic);

% boundary conditions
range = (lds.tps-1)*lds.nphase + (lds.phases);
jac(range,[lds.phases range lds.ncoords+(1:nc)]) = feval(BVP_jac_bc);

