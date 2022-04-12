function [x,v]=Het_adapt_mesh(x,v)
%
% re-interpolates a given x and v on a new mesh such that the
% error in x is approximately uniformly spread over the mesh
% points.
%
global hetds

% construct ups and vps
ups = reshape(x(hetds.coords),hetds.nphase,hetds.tps);
vps = reshape(v(hetds.coords),hetds.nphase,hetds.tps);

% calculate new mesh
tmnew = newmeshHom(ups,hetds.msh,hetds.ntst,hetds.ncol,hetds.ntst,hetds.ncol);

% interpolate ups and vps on new mesh
ups = interp(hetds.msh,hetds.ncol,ups,tmnew,hetds.ncol);
vps = interp(hetds.msh,hetds.ncol,vps,tmnew,hetds.ncol);

% store new mesh
hetds.msh = tmnew;
hetds.finemsh = [0 reshape(repmat(hetds.msh(hetds.tsts),hetds.ncol,1),1,hetds.ntst*hetds.ncol)+kron(hetds.msh((hetds.tsts+1))-hetds.msh(hetds.tsts),((1/hetds.ncol):(1/hetds.ncol):1))];
hetds.dt = hetds.msh(hetds.tsts+1)-hetds.msh(hetds.tsts);

Het_set_ntst_ncol(hetds.ntst,hetds.ncol,tmnew);

x = [ups(:);x(hetds.ncoords+1:end)];
v = [vps(:);v(hetds.ncoords+1:end)];
