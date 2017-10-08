%% function find r,f,g,fdot,gdot,rt,rdot_t,r_d
%
% Author Kuan Chun Hwang

function [f,g,f_dot,g_dot,r,r_t,r_dot_t] = find_r_t_and_r_dot_t(a,phi,dt,r_t0,r_dot_t0,R_t0,R_dot_t0)

% Initialize global variable
global mu_earth;

r = a*(1 - (1 - (R_t0/a))*cos(phi)) + r_t0'*r_dot_t0*sqrt(a/mu_earth)*sin(phi);
f = 1 - (a*(1-cos(phi)))/R_t0;
g = dt - sqrt((a^3)/mu_earth)*(phi - sin(phi));
f_dot = (-sqrt(mu_earth*a)*sin(phi))/(r*R_t0);
g_dot = 1 - ((a*(1-cos(phi)))/r);

% find r(t) and r.(t)
r_t = f*r_t0 + g*r_dot_t0;
r_dot_t = f_dot*r_t0 + g_dot*r_dot_t0;

end