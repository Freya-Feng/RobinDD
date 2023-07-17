function [gx,gy]= specific_g_rb(x,y,g,dy,miny)
%unit_rhs   unit RHS forcing function
%   f = specific_rhs(x,y,nel)
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          nel        number of elements  
%   IFISS function: DJS; 28 February 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
nbc=length(g)/2;
g1=g(1:nbc);g2=g(nbc+1:end);
n=ceil((y-miny)/dy)+1;
gx=g1(n-1)+(y-miny-(n-2)*dy).*(g1(n)-g1(n-1))/dy;
gy=g2(n-1)+(y-miny-(n-2)*dy).*(g2(n)-g2(n-1))/dy;
if max(y)>1
    fprintf('error')
end
return