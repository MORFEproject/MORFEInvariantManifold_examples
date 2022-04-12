function [x,v]=HTHet_new_mesh(x,v,ntst,ncol)
global HTHetds

% some values we need
%oldncoord = HTHetds.ncoords;

% construct ups
ups = reshape(x,HTHetds.nphase,ntst*ncol+1);

% calculate new mesh
tmnew = newmeshHom(ups,HTHetds.msh,ntst,ncol,HTHetds.ntst,HTHetds.ncol);

% interpolate ups on new mesh
ups = interp(HTHetds.msh,ncol,ups,tmnew,HTHetds.ncol);
x = ups(:);

% interpolate vps if supplied
if ~isempty(v)
  vps = interp(HTHetds.msh,ncol,reshape(v,HTHetds.nphase,ntst*ncol+1),tmnew,HTHetds.ncol);
  v = vps(:);
  v = v/norm(v);
end

% store new mesh
HTHet_set_ntst_ncol(HTHetds.ntst,HTHetds.ncol,tmnew);
