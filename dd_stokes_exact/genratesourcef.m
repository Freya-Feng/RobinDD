syms x y
ux=(sin(pi*x))^2*sin(2*pi*y);
uy=-sin(2*pi*x)*(sin(pi*y))^2;

p=cos(pi*x)*sin(pi*y);
fx=-(diff(ux,x,2)+diff(ux,y,2))+diff(p,x)
fy=-(diff(uy,x,2)+diff(uy,y,2))+diff(p,y)