function [x,p] = eqrearr(x0)
global eds
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x) and parameters (p)

% WM: removed loop

l = length(eds.ActiveParams);
ncoo = length(x0)-length(eds.ActiveParams);

p = eds.P0;
p(eds.ActiveParams) = x0(ncoo+(1:l));
x = x0(1:ncoo);
