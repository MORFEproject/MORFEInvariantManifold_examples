function [x,v]=HTHet_adapt_mesh(x,v)
%
% re-interpolates a given x and v on a new mesh such that the
% error in x is approximately uniformly spread over the mesh
% points.
%
global HTHetds

% construct ups and vps
ups = reshape(x(HTHetds.coords),HTHetds.nphase,HTHetds.tps);
vps = reshape(v(HTHetds.coords),HTHetds.nphase,HTHetds.tps);

% calculate new mesh
tmnew = newmeshHom(ups,HTHetds.msh,HTHetds.ntst,HTHetds.ncol,HTHetds.ntst,HTHetds.ncol);

% interpolate ups and vps on new mesh
ups = interp(HTHetds.msh,HTHetds.ncol,ups,tmnew,HTHetds.ncol);
vps = interp(HTHetds.msh,HTHetds.ncol,vps,tmnew,HTHetds.ncol);

% store new mesh
HTHetds.msh = tmnew;
HTHetds.finemsh = [0 reshape(repmat(HTHetds.msh(HTHetds.tsts),HTHetds.ncol,1),1,HTHetds.ntst*HTHetds.ncol)+kron(HTHetds.msh((HTHetds.tsts+1))-HTHetds.msh(HTHetds.tsts),((1/HTHetds.ncol):(1/HTHetds.ncol):1))];
HTHetds.dt = HTHetds.msh(HTHetds.tsts+1)-HTHetds.msh(HTHetds.tsts);

HTHet_set_ntst_ncol(HTHetds.ntst,HTHetds.ncol,tmnew);

x = [ups(:);x(HTHetds.ncoords+1:end)];
v = [vps(:);v(HTHetds.ncoords+1:end)];
