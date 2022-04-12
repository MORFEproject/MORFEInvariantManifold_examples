function w = nc_weight ( n )
%
%  function w = nc_weight ( n )
%
%  For 2 <= N <= 7, returns the weights W of an
%  N+1 point Newton-Cotes quadrature rule.

switch( n )
  case 2
    w = [1/6; 4/6; 1/6];
  case 3
    w = [1/8; 3/8; 3/8; 1/8];
  case 4
    w = [7/90; 32/90; 12/90; 32/90; 7/90];
  case 5
    w = [19/288; 25/96; 25/144; 25/144; 25/96; 19/288];
  case 6
    w = [41/840; 9/35; 9/280; 34/105; 9/280; 9/35; 41/840];
  case 7
    w = [751/17280; 3577/17280; 49/640; 2989/17280; 2989/17280; 49/640; 3577/17280; 751/17280];
end
