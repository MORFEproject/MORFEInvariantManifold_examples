function [x,v]=adapt_mesh2(x,v)
%
% re-interpolates a given x and v on a new mesh such that the
% error in x is approximately uniformly spread over the mesh
% points.
% special version for perioddoubling2
%
global lds

% construct ups and vps

ups1 = reshape(x(1:lds.ncoords),lds.nphase,lds.tps);
vps1 = reshape(v(1:lds.ncoords),lds.nphase,lds.tps);
ups2 = reshape(x(lds.ncoords+1:2*lds.ncoords),lds.nphase,lds.tps);
vps2 = reshape(v(lds.ncoords+1:2*lds.ncoords),lds.nphase,lds.tps);

% calculate new mesh
tmnew = newmesh(ups1,lds.msh,lds.ntst,lds.ncol,lds.ntst,lds.ncol);

% interpolate ups and vps on new mesh
ups1 = interp(lds.msh,lds.ncol,ups1,tmnew,lds.ncol);
ups2 = interp(lds.msh,lds.ncol,ups2,tmnew,lds.ncol);
vps1 = interp(lds.msh,lds.ncol,vps1,tmnew,lds.ncol);
vps2 = interp(lds.msh,lds.ncol,vps2,tmnew,lds.ncol);

% store new mesh
lds.msh = tmnew;
lds.finemsh = [0 reshape(repmat(lds.msh(lds.tsts),lds.ncol,1),1,lds.ntst*lds.ncol)+kron(lds.msh(lds.tsts+1)-lds.msh(lds.tsts),((1/lds.ncol):(1/lds.ncol):1))];
lds.dt = lds.msh(lds.tsts+1)-lds.msh(lds.tsts);

% set result
x = [ups1(:);ups2(:);x(lds.pars)];
v = [vps1(:);vps2(:);v(lds.pars)];
v = v/norm(v); % renormalize v (probably should not be necessery)
