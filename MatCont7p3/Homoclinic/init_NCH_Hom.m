function [x,v] = init_NCH_Hom(odefile, x, v, s, p, ap, ntst, ncol,extravec,T,eps0,eps1)

[x,v] = init_HSN_Hom(odefile, x, v, s, p, ap, ntst, ncol,extravec,T,eps0,eps1);
