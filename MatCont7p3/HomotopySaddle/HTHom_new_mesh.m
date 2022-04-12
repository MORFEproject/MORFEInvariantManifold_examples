function [x,v]=HTHom_new_mesh(x,v,ntst,ncol)
global HTHomds

% some values we need
%oldncoord = HTHomds.ncoords;

% construct ups
ups = reshape(x,HTHomds.nphase,ntst*ncol+1);
% calculate new mesh
tmnew = newmeshHom(ups,HTHomds.msh,ntst,ncol,HTHomds.ntst,HTHomds.ncol);

% interpolate ups on new mesh
ups = interp(HTHomds.msh,ncol,ups,tmnew,HTHomds.ncol);
x = ups(:);

% interpolate vps if supplied
if ~isempty(v)
  vps = interp(HTHomds.msh,ncol,reshape(v,HTHomds.nphase,ntst*ncol+1),tmnew,HTHomds.ncol);
  v = vps(:);
  v = v/norm(v);
end

% store new mesh
HTHom_set_ntst_ncol(HTHomds.ntst,HTHomds.ncol,tmnew);
