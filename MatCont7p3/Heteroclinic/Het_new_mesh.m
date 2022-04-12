function [x,v]=Het_new_mesh(x,v,ntst,ncol)
global hetds

% some values we need
oldncoord = hetds.ncoords;

% construct ups
ups = reshape(x(hetds.coords),hetds.nphase,hetds.tps);

% calculate new mesh
tmnew = newmeshHom(ups,hetds.msh,hetds.ntst,hetds.ncol,ntst,ncol);

% interpolate ups on new mesh
ups = interp(hetds.msh,hetds.ncol,ups,tmnew,ncol);
x = [ups(:);x((oldncoord+1):end)];

% interpolate vps if supplied
if ~isempty(v)
  vps = interp(hetds.msh,hetds.ncol,reshape(v(hetds.coords),hetds.nphase,hetds.tps),tmnew,ncol);
  v = [vps(:);v((oldncoord+1):end)];
  v = v/norm(v);
end

% store new mesh
Het_set_ntst_ncol(ntst,ncol,tmnew);
