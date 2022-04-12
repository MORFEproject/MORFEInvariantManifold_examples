% [T11h, T12h, E21, T22h] = riccatiCoeff(hom, Q0, A, Nsub)
%
% Given Q0 and A, set up the coefficient matrices T22h, T11h, E21, T12h
% for the Riccati equation: T22h*Y - Y*T11h = -E21 + Y*T12h*Y
%
% Inputs:
%  Q0:        an n-by-n (in dense case) old block Schur
%       (with Nsub-dimensional (unstable or stable) invariant subspace)
%  A:         an n-by-n


function [T11h, T12h, E21, T22h] = ricattiCoeff(Q0, A, Nsub)

Th                   = Q0'*A*Q0;

T11h  = Th(1:Nsub,     1:Nsub);
T12h  = Th(1:Nsub,     Nsub+1:end);
E21   = Th(Nsub+1:end, 1:Nsub);
T22h  = Th(Nsub+1:end, Nsub+1:end);
