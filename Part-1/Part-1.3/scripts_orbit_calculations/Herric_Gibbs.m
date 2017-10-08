%% function Herrick- Gibbs
%
% Author: Kuan-Chun Hwang

function orb_parameters = Herric_Gibbs(obs_n1,obs_n2,obs_n3,t_stamp,sat_LGpolar,gs_llh)

% Obtain global constants
global mu_earth;

% Find ground station ECEF
gs_ecef = llhgd2ecef(gs_llh);

% Convert ground station coordinates from ECEF to ECI
gs1_eci = ecef2eci(gs_ecef,obs_n1);
gs2_eci = ecef2eci(gs_ecef,obs_n2);
gs3_eci = ecef2eci(gs_ecef,obs_n3);

% Extract observation data
obs1_LG_polar = sat_LGpolar(:,obs_n1/t_stamp);
obs2_LG_polar = sat_LGpolar(:,obs_n2/t_stamp);
obs3_LG_polar = sat_LGpolar(:,obs_n3/t_stamp);

% Convert from polar to cartesian
obs1_LG = polar2cartesian(obs1_LG_polar);
obs2_LG = polar2cartesian(obs2_LG_polar);
obs3_LG = polar2cartesian(obs3_LG_polar);

% Convert observations from LGDV to ECEF
obs1_ecef = lg2ecef(obs1_LG,gs_llh);
obs2_ecef = lg2ecef(obs2_LG,gs_llh);
obs3_ecef = lg2ecef(obs3_LG,gs_llh);

% Convert relative position of observations from ECEF to ECI
obs1r_eci = ecef2eci(obs1_ecef,obs_n1);
obs2r_eci = ecef2eci(obs2_ecef,obs_n2);
obs3r_eci = ecef2eci(obs3_ecef,obs_n3);

% Find satellite position in ECI
obs1_eci = obs1r_eci + gs1_eci;
obs2_eci = obs2r_eci + gs2_eci;
obs3_eci = obs3r_eci + gs3_eci;

% Find magnitude of each vector
r1_magnitude = norm(obs1_eci);
r2_magnitude = norm(obs2_eci);
r3_magnitude = norm(obs3_eci);

% Define time difference
deltaT32 = (obs_n3 - obs_n2);
deltaT21 = (obs_n2 - obs_n1);
deltaT31 = (obs_n3 - obs_n1);

% Calculate v2
A = mu_earth / (12*(r1_magnitude)^3);
B = mu_earth / (12*(r2_magnitude)^3);
C = mu_earth / (12*(r3_magnitude)^3);

v2 = -deltaT32*((1/(deltaT21*deltaT31))+ A)*obs1_eci + (deltaT32 - deltaT21)*((1/(deltaT21*deltaT32))+B)*obs2_eci + deltaT21*((1/(deltaT32*deltaT31))+C)*obs3_eci;

% Find v_eci
%v_eci_s = v2 + cross(w_earth_eci,gs2_eci);
% Magnitude of velocity vector
v_magnitude = norm(v2);

% Semi-major axis va Vis-Viva Law
%a = ((2/r2_magnitude) - ((v_magnitude^2)/mu_earth))^(-1);
a = (r2_magnitude*mu_earth)/(2*mu_earth - v_magnitude^2*r2_magnitude);
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
capital_omega = atan2(IcrossNdotK,IdotN);

% Find Arugument of perigee
cos_omega = dot(N,e)/norm(e);
sin_omega = dot((cross(N,e)/norm(e)),w);
omega = atan2(sin_omega,cos_omega);

% Find True Anomaly
cos_theta = dot(e,obs2_eci)/(norm(e)*r2_magnitude);
sin_theta = dot(cross(e,obs2_eci)/(norm(e)*r2_magnitude),w);
theta = atan2(sin_theta,cos_theta);

orb_parameters = [a;inclination;e_magni;capital_omega;omega;theta];

end