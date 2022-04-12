function [x,v]=HTHSN_adapt_mesh(x,v)
%
% re-interpolates a given x and v on a new mesh such that the
% error in x is approximately uniformly spread over the mesh
% points.
%
global HTHSNds

% construct ups and vps
ups = reshape(x(HTHSNds.coords),HTHSNds.nphase,HTHSNds.tps);
vps = reshape(v(HTHSNds.coords),HTHSNds.nphase,HTHSNds.tps);

% calculate new mesh
tmnew = newmeshHom(ups,HTHSNds.msh,HTHSNds.ntst,HTHSNds.ncol,HTHSNds.ntst,HTHSNds.ncol);

% interpolate ups and vps on new mesh
ups = interp(HTHSNds.msh,HTHSNds.ncol,ups,tmnew,HTHSNds.ncol);
vps = interp(HTHSNds.msh,HTHSNds.ncol,vps,tmnew,HTHSNds.ncol);

% store new mesh
HTHSNds.msh = tmnew;
HTHSNds.finemsh = [0 reshape(repmat(HTHSNds.msh(HTHSNds.tsts),HTHSNds.ncol,1),1,HTHSNds.ntst*HTHSNds.ncol)+kron(HTHSNds.msh((HTHSNds.tsts+1))-HTHSNds.msh(HTHSNds.tsts),((1/HTHSNds.ncol):(1/HTHSNds.ncol):1))];
HTHSNds.dt = HTHSNds.msh(HTHSNds.tsts+1)-HTHSNds.msh(HTHSNds.tsts);

HTHSN_set_ntst_ncol(HTHSNds.ntst,HTHSNds.ncol,tmnew);

x = [ups(:);x(HTHSNds.ncoords+1:end)];
v = [vps(:);v(HTHSNds.ncoords+1:end)];
