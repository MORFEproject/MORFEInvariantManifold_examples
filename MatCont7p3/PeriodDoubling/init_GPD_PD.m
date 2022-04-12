function [x,v] = init_GPD_PD(odefile, x, s, ap, ntst, ncol)
%
% [x0,v0] = init_GPD_PD(odefile, x, s, ap, ntst, ncol)
%
[x,v] = init_PD_PD(odefile, x, s, ap, ntst, ncol);