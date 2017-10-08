%% function Herrick- Gibbs
%
% Author: Kuan-Chun Hwang

function [orb_parameters, obs2_eci, v2, time_use2, times_obs] = Herric_Gibbs(sat_lg_polar_recieved,gs_llh, satChosen,gs_number,visible_time,times)

% Obtain global constants
global mu_earth;
global w_earth;

% Find ground station ECEF
gs_ecef = llhgc2ecef(gs_llh);

% For each satellite, use the the signal recieved by each ground station to
% find initial orbital parameters
for n = satChosen
    for gs_n = gs_number
        
        % Extract the chosen ground station location
        gs_ecef_use = gs_ecef(:,gs_n);
        
        % Extract the LG position of the chosen satellite to one of the
        % ground station
        sat_lg_vis = sat_lg_polar_recieved;
        
        % Extract the visible times index of the satellite to the ground station
        visTime = visible_time{gs_n,n};
    end
end

% Find the times and take the observation at the index 3,4,5
times_obs = times(visTime);
time_use1 = times_obs(3);
time_use2 = times_obs(4);
time_use3 = times_obs(5);

% Convert ground station coordinates from ECEF to ECI
gs1_eci = ecef2eci(gs_ecef_use,time_use1);
gs2_eci = ecef2eci(gs_ecef_use,time_use2);
gs3_eci = ecef2eci(gs_ecef_use,time_use3);

% Extract observation data
obs1_LG_polar = sat_lg_vis(:,3);
obs2_LG_polar = sat_lg_vis(:,4);
obs3_LG_polar = sat_lg_vis(:,5);

% Convert from polar to cartesian
obs1_LG = polar2cartesian(obs1_LG_polar);
obs2_LG = polar2cartesian(obs2_LG_polar);
obs3_LG = polar2cartesian(obs3_LG_polar);

% Convert observations from LGDV to ECEF
obs1_ecef = lg2ecef(obs1_LG,gs_llh(:,gs_n));
obs2_ecef = lg2ecef(obs2_LG,gs_llh(:,gs_n));
obs3_ecef = lg2ecef(obs3_LG,gs_llh(:,gs_n));

% Convert relative position of observations from ECEF to ECI
obs1r_eci = ecef2eci(obs1_ecef,time_use1);
obs2r_eci = ecef2eci(obs2_ecef,time_use2);
obs3r_eci = ecef2eci(obs3_ecef,time_use3);

% Find satellite position in ECI
obs1_eci = obs1r_eci + gs1_eci;
obs2_eci = obs2r_eci + gs2_eci;
obs3_eci = obs3r_eci + gs3_eci;

% Find magnitude of each vector
r1_magnitude = norm(obs1_eci);
r2_magnitude = norm(obs2_eci);
r3_magnitude = norm(obs3_eci);

% Define time difference
deltaT32 = (time_use3 - time_use2);
deltaT21 = (time_use2 - time_use1);
deltaT31 = (time_use3 - time_use1);

% Calculate v2
A = mu_earth / (12*(r1_magnitude)^3);
B = mu_earth / (12*(r2_magnitude)^3);
C = mu_earth / (12*(r3_magnitude)^3);

v2 = -deltaT32*((1/(deltaT21*deltaT31))+ A)*obs1_eci + (deltaT32 - deltaT21)*((1/(deltaT21*deltaT32))+B)*obs2_eci + deltaT21*((1/(deltaT32*deltaT31))+C)*obs3_eci;

% Find v_eci_s
w_eci_e=[0;0;w_earth];
v_eci = v2 + cross(w_eci_e,gs2_eci);

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