function [x,v]=new_mesh(x,v,ntst,ncol)
global lds

% some values we need
oldncoord = lds.ncoords;

% construct ups
ups = reshape(x(lds.coords),lds.nphase,lds.tps);

% calculate new mesh
tmnew = newmesh(ups,lds.msh,lds.ntst,lds.ncol,ntst,ncol);

% interpolate ups on new mesh
ups = interp(lds.msh,lds.ncol,ups,tmnew,ncol);
x = [ups(:);x((oldncoord+1):end)];

% interpolate vps if supplied
if ~isempty(v)
  vps = interp(lds.msh,lds.ncol,reshape(v(lds.coords),lds.nphase,lds.tps),tmnew,ncol);
  v = [vps(:);v((oldncoord+1):end)];
  v = v/norm(v);
end

% store new mesh
set_ntst_ncol(ntst,ncol,tmnew);
