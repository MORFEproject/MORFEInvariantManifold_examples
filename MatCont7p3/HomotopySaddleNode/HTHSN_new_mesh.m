function [x,v]=HTHSN_new_mesh(x,v,ntst,ncol)
global HTHSNds

% some values we need
%oldncoord = HTHSNds.ncoords;

% construct ups
ups = reshape(x,HTHSNds.nphase,ntst*ncol+1);

% calculate new mesh
tmnew = newmeshHom(ups,HTHSNds.msh,ntst,ncol,HTHSNds.ntst,HTHSNds.ncol);

% interpolate ups on new mesh
ups = interp(HTHSNds.msh,ncol,ups,tmnew,HTHSNds.ncol);
x = ups(:);

% interpolate vps if supplied
if ~isempty(v)
  vps = interp(HTHSNds.msh,ncol,reshape(v,HTHSNds.nphase,ntst*ncol+1),tmnew,HTHSNds.ncol);
  v = vps(:);
  v = v/norm(v);
end

% store new mesh
HTHSN_set_ntst_ncol(HTHSNds.ntst,HTHSNds.ncol,tmnew);
