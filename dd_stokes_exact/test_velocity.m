function [u,v]=test_velocity(x,y)
u=(sin(pi*x)).^2.*sin(2*pi*y);
v=-sin(2*pi*x).*(sin(pi*y)).^2;
end