function [x0,v0] = init_LPNS_LPC(odefile, x, s, ap, ntst, ncol, varargin)
%
% %
% [x0,v0] = init_LPNS_LPC(odefile, x, s, ap, ntst, ncol)
%
%

[x0,v0] = init_LPC_LPC(odefile, x, s, ap, ntst, ncol, varargin{:});