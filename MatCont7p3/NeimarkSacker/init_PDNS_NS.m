function [x0,v0] = init_PDNS_NS(odefile, x, s, ap, ntst, ncol)
%
% [x0,v0] = init_PDNS_NS(odefile, x, s, ap, ntst, ncol)
%
[x0,v0] = init_NS_NS(odefile, x, s, ap, ntst, ncol);