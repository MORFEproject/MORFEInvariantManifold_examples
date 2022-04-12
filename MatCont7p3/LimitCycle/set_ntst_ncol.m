function set_ntst_ncol(ntst,ncol,newmsh)
%
% This function sets the number of mesh and collocation points
% as well as a new mesh distribution. All things that depend on
% the number of mesh and/or collocation points are updated.
%
global cds lds
lds.ntst = ntst;
lds.ncol = ncol;

lds.tsts = 1:ntst;
lds.cols = 1:ncol;
lds.tps = lds.ntst*lds.ncol+1; % number of points on curve

lds.ncoords = lds.tps*lds.nphase;
lds.coords = 1:lds.ncoords;
lds.PeriodIdx = lds.ncoords+1;
lds.idxmat = reshape(fix((1:((lds.ncol+1)*lds.ntst))/(1+1/lds.ncol))+1,lds.ncol+1,lds.ntst);
cds.ndim = lds.ncoords+2;
cds.oldJacX = [];

lds.msh = newmsh;
lds.finemsh = [0 reshape(repmat(lds.msh(1:lds.ntst),lds.ncol,1),1,lds.ntst*lds.ncol)+kron(lds.msh(2:(lds.ntst+1))-lds.msh(1:lds.ntst),((1/lds.ncol):(1/lds.ncol):1))];
lds.dt = lds.msh(lds.tsts+1)-lds.msh(lds.tsts);

lds.upoldp = [];
lds.multipliersX = [];
lds.CalcMultipliers =0;
calc_weigths;
