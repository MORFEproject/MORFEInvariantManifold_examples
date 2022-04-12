function HTHSN_set_ntst_ncol(ntst,ncol,newmsh)
%
% This function sets the number of mesh and collocation points
% as well as a new mesh distribution. All things that depend on
% the number of mesh and/or collocation points are updated.
%
global cds HTHSNds

HTHSNds.ntst = ntst;
HTHSNds.ncol = ncol;

HTHSNds.tsts = 1:ntst;
HTHSNds.cols = 1:ncol;
HTHSNds.tps = HTHSNds.ntst*HTHSNds.ncol+1; % number of points on curve

HTHSNds.ncoords = HTHSNds.tps*HTHSNds.nphase;
HTHSNds.coords = 1:HTHSNds.ncoords;
HTHSNds.PeriodIdx = HTHSNds.ncoords+HTHSNds.nphase;
HTHSNds.idxmat = reshape(fix((1:((HTHSNds.ncol+1)*HTHSNds.ntst))/(1+1/HTHSNds.ncol))+1,HTHSNds.ncol+1,HTHSNds.ntst);
if isfield(HTHSNds,'ndim')
    cds.ndim = HTHSNds.ndim;
end
cds.oldJacX = [];

HTHSNds.msh = newmsh;
HTHSNds.finemsh = [0 reshape(repmat(HTHSNds.msh(1:HTHSNds.ntst),HTHSNds.ncol,1),1,HTHSNds.ntst*HTHSNds.ncol)+kron(HTHSNds.msh(2:(HTHSNds.ntst+1))-HTHSNds.msh(1:HTHSNds.ntst),((1/HTHSNds.ncol):(1/HTHSNds.ncol):1))];
HTHSNds.dt = HTHSNds.msh(HTHSNds.tsts+1)-HTHSNds.msh(HTHSNds.tsts);

HTHSNds.upold = [];
HTHSNds.upoldp = [];
HTHSN_calc_weigths;
