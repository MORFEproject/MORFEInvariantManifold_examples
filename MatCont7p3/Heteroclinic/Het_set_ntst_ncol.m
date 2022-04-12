function Het_set_ntst_ncol(ntst,ncol,newmsh)
%
% This function sets the number of mesh and collocation points
% as well as a new mesh distribution. All things that depend on
% the number of mesh and/or collocation points are updated.
%
global cds hetds

hetds.ntst = ntst;
hetds.ncol = ncol;

hetds.tsts = 1:ntst;
hetds.cols = 1:ncol;
hetds.tps = hetds.ntst*hetds.ncol+1; % number of points on curve

hetds.ncoords = hetds.tps*hetds.nphase;
hetds.coords = 1:hetds.ncoords;
hetds.PeriodIdx = hetds.ncoords+2*hetds.nphase;
hetds.idxmat = reshape(fix((1:((hetds.ncol+1)*hetds.ntst))/(1+1/hetds.ncol))+1,hetds.ncol+1,hetds.ntst);
if isfield(hetds,'ndim')
    cds.ndim = hetds.ndim;
end
cds.oldJacX = [];

hetds.msh = newmsh;
hetds.finemsh = [0 reshape(repmat(hetds.msh(1:hetds.ntst),hetds.ncol,1),1,hetds.ntst*hetds.ncol)+kron(hetds.msh(2:(hetds.ntst+1))-hetds.msh(1:hetds.ntst),((1/hetds.ncol):(1/hetds.ncol):1))];
hetds.dt = hetds.msh(hetds.tsts+1)-hetds.msh(hetds.tsts);

hetds.upold = [];
hetds.upoldp = [];
Het_calc_weigths;
