function [x0,v0]= init_HH_H(odefile, x, p, ap)
%
% [x1,v1] = init_HH_H(odefile, x, p, ap)
%
% Initializes a Hopf bifurcation continuation from a Hopf point
% 
%

[x0,v0] = init_H_H(odefile, x, p, ap);