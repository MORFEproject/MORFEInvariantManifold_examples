function [x0,v0] = init_R4_NS(odefile, x, s, ap, ntst, ncol)
%
% [x0,v0] = init_R4_NS(odefile, x, s, ap, ntst, ncol)
%
[x0,v0] = init_NS_NS(odefile, x, s, ap, ntst, ncol);