function HTHet_set_ntst_ncol(ntst,ncol,newmsh)
%
% This function sets the number of mesh and collocation points
% as well as a new mesh distribution. All things that depend on
% the number of mesh and/or collocation points are updated.
%
global cds HTHetds

HTHetds.ntst = ntst;
HTHetds.ncol = ncol;

HTHetds.tsts = 1:ntst;
HTHetds.cols = 1:ncol;
HTHetds.tps = HTHetds.ntst*HTHetds.ncol+1; % number of points on curve

HTHetds.ncoords = HTHetds.tps*HTHetds.nphase;
HTHetds.coords = 1:HTHetds.ncoords;
HTHetds.PeriodIdx = HTHetds.ncoords+2*HTHetds.nphase;
HTHetds.idxmat = reshape(fix((1:((HTHetds.ncol+1)*HTHetds.ntst))/(1+1/HTHetds.ncol))+1,HTHetds.ncol+1,HTHetds.ntst);
if isfield(HTHetds,'ndim')
    cds.ndim = HTHetds.ndim;
end
cds.oldJacX = [];

HTHetds.msh = newmsh;
HTHetds.finemsh = [0 reshape(repmat(HTHetds.msh(1:HTHetds.ntst),HTHetds.ncol,1),1,HTHetds.ntst*HTHetds.ncol)+kron(HTHetds.msh(2:(HTHetds.ntst+1))-HTHetds.msh(1:HTHetds.ntst),((1/HTHetds.ncol):(1/HTHetds.ncol):1))];
HTHetds.dt = HTHetds.msh(HTHetds.tsts+1)-HTHetds.msh(HTHetds.tsts);

HTHetds.upold = [];
HTHetds.upoldp = [];
HTHet_calc_weigths;
