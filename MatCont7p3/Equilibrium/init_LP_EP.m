function [x0,v0]= init_LP_EP(odefile, x, p, ap,varargin)
%
%[x0,v0]= init_LP_EP(odefile, x, p, ap)
%
% Defines odefile
% Sets all parameters for an equilibrium continuation (p)
% and the active parameter (ap)
%
[x0,v0] = init_EP_EP(odefile, x, p, ap);