function [x0,v0]= init_CP_LP(odefile, x, p, ap)
%
% [x0,v0] = init_CP_LP(odefile, x, p, ap) 
% 
%
[x0,v0]= init_LP_LP(odefile, x, p, ap);