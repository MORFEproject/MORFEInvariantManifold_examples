function p = gl_pos ( n )
%
%  function p = gl_pos ( n )
%
%  For 2 <= N <= 7
%  N point Gauss-Legendre quadrature rule over the interval [-1,1].

switch( n )
  case 2
    c = 1/sqrt(3);
    p = [-c; c];
  case 3
    c = sqrt(6/10);
    p = [-c; 0; c];
  case 4
    c = sqrt(24/245);
    c1 = sqrt(3/7+c);
    c2 = sqrt(3/7-c);
    p = [-c1; -c2; c2; c1];
  case 5
    c1 = 0.90617984593866399280;
    c2 = 0.53846931010568309104;
    p = [-c1; -c2; 0; c2; c1];
  case 6
    c1 = 0.93246951420315202781;
    c2 = 0.66120938646626451366;
    c3 = 0.23861918608319690863;
    p = [-c1; -c2; -c3; c3; c2; c1];
  case 7
    c1 = 0.949107991234275852452;
    c2 = 0.74153118559939443986;
    c3 = 0.40584515137739716690;
    p = [-c1; -c2; -c3; 0; c3; c2; c1];
end
