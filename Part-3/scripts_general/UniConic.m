%% Universal Conic Section Method
% This function implements Universal Conic Section Method
%
% Author: Kuan Chun Hwang

function  [ r1,v1 ] = UniConic( r0,v0, dE )

% Initialize global variable
global mu_earth;

% Compute semi-Major axis using vis viva law
a = ((2/norm(r0)) - ((norm(v0)^2)/mu_earth))^-1;

% Calculate dt using the relevant equation provided
dt = sqrt((a^3)/mu_earth)*(dE - (1 - (norm(r0)/a))*sin(dE) + (r0'*v0)/sqrt(mu_earth*a)*(1-cos(dE)));

[f,g,fDot,gDot] = find_f_g_fDot_gDot(a,dE,dt,r0,v0);

% Calculate r1 and v1
r1=f*r0+g*v0;
v1=fDot*r0+gDot*v0;

end