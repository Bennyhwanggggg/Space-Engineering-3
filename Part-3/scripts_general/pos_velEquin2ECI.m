%%  Convert position vector and velocity vector from equinoctial states to ECI
%
% Author: Kuan Chun Hwang

function [pos_vector, vel_vector] = pos_velEquin2ECI(x_parameters)

% Initialize global variable
global mu_earth

% Extract Parameter values
p = x_parameters(1,:);
f = x_parameters(2,:);
g = x_parameters(3,:);
h = x_parameters(4,:);
k = x_parameters(5,:);
L = x_parameters(6,:);

% Find r, s_sq, w, alpha_sq
w = 1 + f.*cos(L) + g.*sin(L);
s_sq = 1 + (h^2) + (k^2);
r = p ./ w;
alpha2 = h^2 - k^2;

% Convert Equinoctial elements to ECI position vector
pos_vector = [ (r/(s_sq))*(cos(L)+alpha2*cos(L)+2*h*k*sin(L));
    (r/(s_sq))*(sin(L)-alpha2*sin(L)+2*h*k*cos(L));
    ((2*r)/(s_sq))*(h*sin(L)-k*cos(L))];

vel_vector = [(-1/s_sq^2)*sqrt(mu_earth/p)*(sin(L)+alpha2*sin(L)-2*h*k*cos(L)+g-2*f*h*k+alpha2*g);
            (-1/s_sq^2)*sqrt(mu_earth/p)*(-cos(L)+alpha2*cos(L)+2*h*k*sin(L)-f+2*g*h*k+alpha2*f);
            (2/(s_sq^2))*sqrt(mu_earth/p)*(h*cos(L)+k*sin(L)+f*h+g*k)];

end