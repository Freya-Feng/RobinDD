function f = specific_rhs_rb(x,y)
%unit_rhs   unit RHS forcing function
%   f = specific_rhs(x,y,nel)
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          nel        number of elements  
%   IFISS function: DJS; 28 February 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
fx=6*pi^2*sin(pi*x).^2.*sin(2*pi*y) - 2*pi^2*cos(pi*x).^2.*sin(2*pi*y) - pi*sin(pi*x).*sin(pi*y);
fy=2*pi^2*cos(pi*y).^2.*sin(2*pi*x) - 6*pi^2*sin(2*pi*x).*sin(pi*y).^2 + pi*cos(pi*x).*cos(pi*y);
f=[fx;fy];
return