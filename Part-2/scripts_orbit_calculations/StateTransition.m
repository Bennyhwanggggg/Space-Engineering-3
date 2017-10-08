%% function Compute state transition Jacobian matrix
%
% Author Kuan Chun Hwang

function dxt_dxt0 = StateTransition(dt, a, f, g, f_dot, g_dot, r_t0, r_dot_t0, R_t0, R_dot_t0, r_t, r_dot_t, R_t, R_dot_t)

% Initialize global variable
global mu_earth;

% Compute relevant variables
alpha = 1/a;
chi = alpha*sqrt(mu_earth)*dt + (r_t'*r_dot_t - r_t0'*r_dot_t0) / sqrt(mu_earth);

u2 = (1 - cos(sqrt(alpha)*chi)) / alpha;
u3 = (sqrt(alpha)*chi - sin(sqrt(alpha)*chi)) / (alpha*sqrt(alpha));
u4 = (chi^2)/(2*alpha) - u2/alpha;
u5 = (chi^3)/(6*alpha) - u3/alpha;

c = (3*u5 - chi*u4 - sqrt(mu_earth)*dt*u2) / sqrt(mu_earth);

% Compute the matrix components
p11 = (R_t/mu_earth)*(r_dot_t - r_dot_t0)*((r_dot_t - r_dot_t0)') + (R_t0^-3)*(R_t0*(1-f)*r_t*r_t0'+c*r_dot_t*r_t0') + f*eye(3);
p12 = (R_t0/mu_earth)*(1-f)*((r_t-r_t0)*r_dot_t0' - (r_dot_t - r_dot_t0)*r_t0') + (c/mu_earth)*r_dot_t*r_dot_t0' + g*eye(3);
p21 = -(R_t0^-2)*(r_dot_t-r_dot_t0)*r_t0' - (R_t^-2)*r_t*(r_dot_t-r_dot_t0)' - ((mu_earth*c)/((R_t^3)*(R_t0^3)))*r_t*r_t0' + f_dot*(eye(3) - (R_t0^-2)*r_t*r_t' + (1/(mu_earth*R_t))*(r_t*r_dot_t' - r_dot_t*r_t')*r_t*(r_dot_t - r_dot_t0)');
p22 = (R_t0/mu_earth)*(r_dot_t - r_dot_t0)*(r_dot_t - r_dot_t0)' + (R_t0^-3)*(R_t*(1-f)*r_t*r_t0' - c*r_t*r_dot_t0') + g_dot*eye(3);

dxt_dxt0 = [p11 p12; p21 p22];
end