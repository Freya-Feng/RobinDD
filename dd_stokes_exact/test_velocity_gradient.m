function [dudx,dudy,dvdx,dvdy]=test_velocity_gradient(x,y)
dudx=2*pi*cos(pi*x).*sin(pi*x).*sin(2*pi*y);
dudy=2*pi*cos(2*pi*y).*sin(pi*x).^2;
dvdx=-2*pi.*cos(2*pi*x).*sin(pi*y).^2;
dvdy=-2*pi.*cos(pi*y).*sin(2*pi*x).*sin(pi*y);
end