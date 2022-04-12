function Hom_set_ntst_ncol(ntst,ncol,nwmesh)
%
% This function sets the number of mesh and collocation points
% as well as a new mesh distribution. All things that depend on
% the number of mesh and/or collocation points are updated.
%
global cds homds
homds.ntst = ntst;
homds.ncol = ncol;

homds.tsts = 1:ntst;
homds.cols = 1:ncol;
homds.tps = homds.ntst*homds.ncol+1; % number of points on curve

homds.ncoords = homds.tps*homds.nphase;
homds.coords = 1:homds.ncoords;
homds.PeriodIdx = homds.ncoords+homds.nphase;
homds.idxmat = reshape(fix((1:((homds.ncol+1)*homds.ntst))/(1+1/homds.ncol))+1,homds.ncol+1,homds.ntst);
if isfield(homds,'ndim')
    cds.ndim = homds.ndim;
end
cds.oldJacX = [];

homds.msh = nwmesh;
homds.finemsh = [0 reshape(repmat(homds.msh(1:homds.ntst),homds.ncol,1),1,homds.ntst*homds.ncol)+kron(homds.msh(2:(homds.ntst+1))-homds.msh(1:homds.ntst),((1/homds.ncol):(1/homds.ncol):1))];
homds.dt = homds.msh(homds.tsts+1)-homds.msh(homds.tsts);

homds.upold = [];
homds.upoldp = [];
Hom_calc_weights;
