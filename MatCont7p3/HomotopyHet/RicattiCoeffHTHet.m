% [T11h, T12h, E21, T22h] = RiccatiCoeff( Q0, A, Nsub)
%
% Given Q0 and A, set up the coefficient matrices R11, R12, E21, R22
% for the Riccati equation: R22*Y - Y*R11 = -E21 + Y*R12*Y
%
% Inputs:
%  Q0:        an n-by-n (in dense case) old block Schur
%       (with Nsub-dimensional (unstable or stable) invariant subspace)
%  A:         an n-by-n Jacobian matrix of the fixed point x0 or x1


function [T11h, T12h, E21, T22h] = RicattiCoeffHTHet(Q0, A, Nsub)

Th                   = Q0'*A*Q0;
T11h  = Th(1:Nsub,     1:Nsub);
T12h  = Th(1:Nsub,     Nsub+1:end);
E21   = Th(Nsub+1:end, 1:Nsub);
T22h  = Th(Nsub+1:end, Nsub+1:end);
