%% Fuel Burn Dynamic Model
% This function runs through the dynamic model by using Universal Conic
% Section Method to obtain r and rDot vector of initial cost, deltaV of
% first burn, coast into final Orbit, entering final orbit and forms the
% relevant equations: objective function, constraint function and
% the augmented lagrangian
%
% Author: Kuan Chun Hwang

function [L,c,f] = BurnDynamicModel(x_parameter,r0,v0,scaleFactor,lambda,rho)

% Obtain global constant
global mu_earth;
global v_initial;
global v_final;
global a_i;
global a_f;
global i_f_rad;
global RAAN_f_rad;
global rDot_final;
global r_final;
% Extract varaibles
x_parameter = x_parameter./scaleFactor;
dE1 = x_parameter(1);
dV1 = x_parameter(2);
theta1 = x_parameter(3);
psi1 = x_parameter(4);
dE2 = x_parameter(5);
dV2 = x_parameter(6);
theta2 = x_parameter(7);
psi2 = x_parameter(8);

% Satellite coasting in initial parking orbit before burn
% Find the vector r1 and v1 which is the position and velocity vector after
% the first change in eccentric anomaly (delta E = pi)
[ r1,v1 ] = UniConic(r0, v0, dE1);

% First Burn
% Calculate the change in velocity required for first burn
dV1_ECI = deltaV(dV1,v1,r1,theta1,psi1); 

% Compute v2
v2 = v1 + dV1_ECI;

% Step3: Moving towards GEO orbit
% Assume instantaneous burn
r2 = r1;

% Repeat universal conic method
[ r3,v3 ] = UniConic( r2,v2, dE2);

% Transit into GEO orbit
% Assume instantaneous burn
r4 = r3; 

% Calculate the change in velocity required for second burn to transit into
% desired orbit
dV2_ECI = deltaV(dV2,v3,r4,theta2,psi2);
v4 = v3 + dV2_ECI;

% Compute the output functions
% Constraint Function
c = [norm(r4)/a_f - 1; norm(v4)/v_final - 1; (v4'*r4)/norm(v4)/norm(r4); v4(3)/norm(v4)-rDot_final(3)/norm(rDot_final); r4(3)/norm(r4)-r_final(3)/norm(r_final)];

% Objective Function deltav1 + deltav2
f = scaleFactor(2)*norm(dV1_ECI) + scaleFactor(6)*norm(dV2_ECI);   

% Augmented lagrangian
L = f - lambda'*c + (rho)*(c')*c;
end