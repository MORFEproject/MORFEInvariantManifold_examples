function [x,v] = init_NCH_HSN(odefile, x, v, s, p, ap, ntst, ncol,extravec,T,eps0,eps1)

[x,v] = init_Hom_HSN(odefile, x, v, s, p, ap, ntst, ncol,extravec,T,eps0,eps1);