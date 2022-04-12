function [x,v]=adapt_mesh(x,v)
%
% re-interpolates a given x and v on a new mesh such that the
% error in x is approximately uniformly spread over the mesh
% points.
%
global lds

% construct ups and vps
ups = reshape(x(lds.coords),lds.nphase,lds.tps);
vps = reshape(v(lds.coords),lds.nphase,lds.tps);

% calculate new mesh
tmnew = newmesh(ups,lds.msh,lds.ntst,lds.ncol,lds.ntst,lds.ncol);

% interpolate ups and vps on new mesh
ups = interp(lds.msh,lds.ncol,ups,tmnew,lds.ncol);
vps = interp(lds.msh,lds.ncol,vps,tmnew,lds.ncol);

% store new mesh
lds.msh = tmnew;
%lds.finemsh = [0 reshape(repmat(lds.msh(lds.tsts),lds.ncol,1),1,lds.ntst*lds.ncol)+kron(lds.msh(lds.tsts)-lds.msh(lds.tsts),((1/lds.ncol):(1/lds.ncol):1))];
lds.finemsh = [0 reshape(repmat(lds.msh(lds.tsts),lds.ncol,1),1,lds.ntst*lds.ncol)+kron(lds.msh((lds.tsts+1))-lds.msh(lds.tsts),((1/lds.ncol):(1/lds.ncol):1))];

lds.dt = lds.msh(lds.tsts+1)-lds.msh(lds.tsts);

% set result
x = [ups(:);x(lds.pars)];
v = [vps(:);v(lds.pars)];
v = v/norm(v); % renormalize v (probably should not be necessery)
