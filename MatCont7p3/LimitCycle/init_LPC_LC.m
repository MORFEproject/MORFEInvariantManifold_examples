function [x,v] = init_LPC_LC(odefile, x, v, s, par, ap, ntst, ncol)
%
% [x0,v0] = init_LPC_LC(odefile, x, v, s, ap, ntst, ncol)
%
% Initializes a limit cycle continuation from a cycle calculated
% in a previous run.
%
global lds

% Make sure size of x and v vectors is right for calling init_LC_LC
x = x(1:lds.ncoords+2,:);
v = v(1:lds.ncoords+2,:);
[x,v] = init_LC_LC(odefile, x, v, s, par, ap, ntst, ncol);