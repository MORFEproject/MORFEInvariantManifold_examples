function [x0,v0] = init_CPC_LPC(odefile, x, s, ap, ntst, ncol, varargin)
%
% %
% [x0,v0] = init_CPC_LPC(odefile, x, s, ap, ntst, ncol)
%
%

[x0,v0] = init_LPC_LPC(odefile, x, s, ap, ntst, ncol, varargin{:});