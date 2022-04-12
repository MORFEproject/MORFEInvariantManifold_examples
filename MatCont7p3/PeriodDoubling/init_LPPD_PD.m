function [x,v] = init_LPPD_PD(odefile, x, s, ap, ntst, ncol)
%
% [x0,v0] = init_LPPD_PD(odefile, x, s, ap, ntst, ncol)
%
[x,v] = init_PD_PD(odefile, x, s, ap, ntst, ncol);