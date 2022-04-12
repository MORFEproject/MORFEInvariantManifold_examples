function [x,v]=Hom_new_mesh(x,v,ntst,ncol)
global homds

% some values we need
oldncoord = homds.ncoords;

% construct ups
ups = reshape(x(homds.coords),homds.nphase,homds.tps);

% calculate new mesh
tmnew = newmeshHom(ups,homds.msh,homds.ntst,homds.ncol,ntst,ncol);

% interpolate ups on new mesh
ups = interp(homds.msh,homds.ncol,ups,tmnew,ncol);
x = [ups(:);x((oldncoord+1):end)];

% interpolate vps if supplied
if ~isempty(v)
  vps = interp(homds.msh,homds.ncol,reshape(v(homds.coords),homds.nphase,homds.tps),tmnew,ncol);
  v = [vps(:);v((oldncoord+1):end)];
  v = v/norm(v);
end

% store new mesh
Hom_set_ntst_ncol(ntst,ncol,tmnew);

