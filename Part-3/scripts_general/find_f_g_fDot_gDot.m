%%  Find f,g,fdot,gdot
%
% Author: Kuan Chun Hwang

function [f,g,fDot,gDot] = find_f_g_fDot_gDot(a,phi,dt,r0,v0)

% Initialize global variable
global mu_earth;

r = a*(1 - (1 - (norm(r0)/a))*cos(phi)) + r0'*v0*sqrt(a/mu_earth)*sin(phi);
f = 1 - (a*(1-cos(phi)))/norm(r0);
g = dt - sqrt((a^3)/mu_earth)*(phi - sin(phi));
fDot = (-sqrt(mu_earth*a)*sin(phi))/(r*norm(r0));
gDot = 1 - ((a*(1-cos(phi)))/r);

end