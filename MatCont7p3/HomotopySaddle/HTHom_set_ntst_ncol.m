function HTHom_set_ntst_ncol(ntst,ncol,newmsh)
%
% This function sets the number of mesh and collocation points
% as well as a new mesh distribution. All things that depend on
% the number of mesh and/or collocation points are updated.
%
global cds HTHomds

HTHomds.ntst = ntst;
HTHomds.ncol = ncol;

HTHomds.tsts = 1:ntst;
HTHomds.cols = 1:ncol;
HTHomds.tps = HTHomds.ntst*HTHomds.ncol+1; % number of points on curve

HTHomds.ncoords = HTHomds.tps*HTHomds.nphase;
HTHomds.coords = 1:HTHomds.ncoords;
HTHomds.PeriodIdx = HTHomds.ncoords+HTHomds.nphase;
HTHomds.idxmat = reshape(fix((1:((HTHomds.ncol+1)*HTHomds.ntst))/(1+1/HTHomds.ncol))+1,HTHomds.ncol+1,HTHomds.ntst);
if isfield(HTHomds,'ndim')
    cds.ndim = HTHomds.ndim;
end
cds.oldJacX = [];

HTHomds.msh = newmsh;
HTHomds.finemsh = [0 reshape(repmat(HTHomds.msh(1:HTHomds.ntst),HTHomds.ncol,1),1,HTHomds.ntst*HTHomds.ncol)+kron(HTHomds.msh(2:(HTHomds.ntst+1))-HTHomds.msh(1:HTHomds.ntst),((1/HTHomds.ncol):(1/HTHomds.ncol):1))];
HTHomds.dt = HTHomds.msh(HTHomds.tsts+1)-HTHomds.msh(HTHomds.tsts);

HTHomds.upold = [];
HTHomds.upoldp = [];
HTHom_calc_weigths;
