function [x,v]=Hom_adapt_mesh(x,v)
%
% re-interpolates a given x and v on a new mesh such that the
% error in x is approximately uniformly spread over the mesh
% points.
%
global homds

% construct ups and vps
ups = reshape(x(homds.coords),homds.nphase,homds.tps);
vps = reshape(v(homds.coords),homds.nphase,homds.tps);

% calculate new mesh
tmnew = newmeshHom(ups,homds.msh,homds.ntst,homds.ncol,homds.ntst,homds.ncol);

% interpolate ups and vps on new mesh
ups = interp(homds.msh,homds.ncol,ups,tmnew,homds.ncol);
vps = interp(homds.msh,homds.ncol,vps,tmnew,homds.ncol);

% store new mesh
homds.msh = tmnew;
homds.finemsh = [0 reshape(repmat(homds.msh(homds.tsts),homds.ncol,1),1,homds.ntst*homds.ncol)+kron(homds.msh((homds.tsts+1))-homds.msh(homds.tsts),((1/homds.ncol):(1/homds.ncol):1))];

homds.dt = homds.msh(homds.tsts+1)-homds.msh(homds.tsts);

Hom_set_ntst_ncol(homds.ntst,homds.ncol,tmnew);

x = [ups(:);x(homds.ncoords+1:end)];
v = [vps(:);v(homds.ncoords+1:end)];
