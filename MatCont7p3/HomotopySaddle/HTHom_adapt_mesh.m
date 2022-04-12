function [x,v]=HTHom_adapt_mesh(x,v)
%
% re-interpolates a given x and v on a new mesh such that the
% error in x is approximately uniformly spread over the mesh
% points.
%
global HTHomds

% construct ups and vps
ups = reshape(x(HTHomds.coords),HTHomds.nphase,HTHomds.tps);
vps = reshape(v(HTHomds.coords),HTHomds.nphase,HTHomds.tps);

% calculate new mesh
tmnew = newmeshHom(ups,HTHomds.msh,HTHomds.ntst,HTHomds.ncol,HTHomds.ntst,HTHomds.ncol);

% interpolate ups and vps on new mesh
ups = interp(HTHomds.msh,HTHomds.ncol,ups,tmnew,HTHomds.ncol);
vps = interp(HTHomds.msh,HTHomds.ncol,vps,tmnew,HTHomds.ncol);

% store new mesh
HTHomds.msh = tmnew;
HTHomds.finemsh = [0 reshape(repmat(HTHomds.msh(HTHomds.tsts),HTHomds.ncol,1),1,HTHomds.ntst*HTHomds.ncol)+kron(HTHomds.msh((HTHomds.tsts+1))-HTHomds.msh(HTHomds.tsts),((1/HTHomds.ncol):(1/HTHomds.ncol):1))];
HTHomds.dt = HTHomds.msh(HTHomds.tsts+1)-HTHomds.msh(HTHomds.tsts);

HTHom_set_ntst_ncol(HTHomds.ntst,HTHomds.ncol,tmnew);

x = [ups(:);x(HTHomds.ncoords+1:end)];
v = [vps(:);v(HTHomds.ncoords+1:end)];
