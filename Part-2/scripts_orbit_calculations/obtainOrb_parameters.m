%% function Herrick- Gibbs
%
% Author: Kuan-Chun Hwang

function [orb_parameters] = obtainOrb_parameters(x0)

% Obtain global constants
global mu_earth;
global w_earth;

r2 = x0(:,1);
v2 = x0(:,2);

% Retrive satellite position in ECI
obs2_eci = r2;

% Find magnitude of each vector
r2_magnitude = norm(obs2_eci);

% Magnitude of velocity vector
v_magnitude = norm(v2);

% Semi-major axis va Vis-Viva Law with auto correction
a = (r2_magnitude*mu_earth)/(2*mu_earth - (v_magnitude^2)*r2_magnitude);

% Find inclincation
w = cross(obs2_eci,v2)/norm(cross(obs2_eci,v2));
K = [0 0 1];
inclination = acos(dot(w,K));

% find eccentricity
e = (1/mu_earth)*(((v_magnitude^2)-(mu_earth/r2_magnitude))*obs2_eci - dot(obs2_eci,v2)*v2);
e_magni = norm(e);

% Find right ascension of the ascending node
N = cross(K,w)/norm(cross(K,w));
I = [1 0 0];
IdotN = dot(I,N);
IcrossNdotK = dot(cross(I,N),K);
captial_omega = atan2(IcrossNdotK,IdotN);

% Find Arugument of perigee
cos_omega = dot(N,e)/norm(e);
sin_omega = dot((cross(N,e)/norm(e)),w);
omega = atan2(sin_omega,cos_omega);

% Find True Anomaly
cos_theta = dot(e,obs2_eci)/(norm(e)*r2_magnitude);
sin_theta = dot(cross(e,obs2_eci)/(norm(e)*r2_magnitude),w);
theta = atan2(sin_theta,cos_theta);


% Take the average of the parameters found for each ground station and
% use this as the initial parameter
orb_parameters = [a;inclination;e_magni;captial_omega;omega;theta];

end