%% function 'Ground Trace of satellite'
%
% Plots the orbit of satellite in ECEF
%
% Input  :   inclination, right_asc, eccentricity, omega,
%            mean_anomaly, mean_motion, semi_major_axis
% Outputs:  3D plot of orbit in ECEF position
%
% Kuan Chun Hwang

%% This orbital calculation function will plot 3D orbit plot and ground tracking plot
function ground_trace(inclination, right_asc, eccentricity, omega, mean_anomaly, mean_motion,t_stamp)

% Initialize global vairiables
global mu_earth;

% Convert Degrees to radians
inclination_rad = inclination * pi / 180;
right_asc_rad = right_asc * pi / 180;
mean_anomaly_rad = mean_anomaly * pi / 180;
omega_rad = omega * pi / 180;

% Find semi major axis
semi_major_axis = nthroot(mu_earth/(mean_motion^2),3); %(m)

% Initialize first guess of eccentric anomaly which is t0 = t so mean
% anomaly is the same as eccentric anomaly
eccentric_anomaly_old_rad = mean_anomaly_rad;

% Set tolerance which is the accuracy of the iterations and assume
% difference is great at the start
toler = 0.0000000000000001;
delta= 10000000;

position_number = 1;

% Orbit plotting loop

for t=0:t_stamp:86400
    
    % Set mean anomay at any given time in that time step
    mean_anomaly_rad = mean_anomaly_rad + (sqrt(mu_earth/(semi_major_axis^3)))*t_stamp;
    
    % Start Iteration
    while delta > toler
        
        % Solve Keplers equation for E
        delta = eccentric_anomaly_old_rad - eccentricity * sin(eccentric_anomaly_old_rad) - mean_anomaly_rad;
        
        % Next iteration value
        eccentric_anomaly_new_rad = eccentric_anomaly_old_rad - delta / (1 - eccentricity * cos(eccentric_anomaly_old_rad));
        eccentric_anomaly_old_rad = eccentric_anomaly_new_rad;
        
    end
    
    % reset delta value
    delta = 10000000000000;
    
    % Final eccentric anomaly value
    essentric_anomaly_ans = eccentric_anomaly_old_rad;
    
    % Find True Anomaly
    true_anomaly_rad = 2 * atan( sqrt((1+eccentricity)/(1-eccentricity)) * tan (essentric_anomaly_ans/2));
    
    % Find semilatus rectum (p)
    semilatus_rectum = semi_major_axis * (1 - eccentricity^2);
    
    % Calculate Distance to satellite
    distance = semilatus_rectum ./ (1 + eccentricity * cos (true_anomaly_rad));
    
    % Coordinates in orbital frame
    x_orbital = distance * cos(true_anomaly_rad);
    y_orbital = distance * sin(true_anomaly_rad);
    z_orbital = 0;
    
    satellite_perifocal_position = [x_orbital; y_orbital; z_orbital];
    
    % Calculate velocities
    vx_orbital = -sqrt(mu_earth/semilatus_rectum) * sin(true_anomaly_rad);
    vy_orbital = -sqrt(mu_earth/semilatus_rectum) * (eccentricity + cos(true_anomaly_rad));
    vz_orbital = 0;
    
    satellite_perifocal_vel = [vx_orbital; vy_orbital; vz_orbital];
    
    % Initialize Transfer matrix from orbital coordinates to ECI frame
    orbit_to_ECI = [cos(right_asc_rad)*cos(omega_rad)-sin(right_asc_rad)*sin(omega_rad)*cos(inclination_rad) -cos(right_asc_rad)*sin(omega_rad)-sin(right_asc_rad)*cos(omega_rad)*cos(inclination_rad) sin(right_asc_rad)*sin(inclination_rad)
        sin(right_asc_rad)*cos(omega_rad)+cos(right_asc_rad)*sin(omega_rad)*cos(inclination_rad) -sin(right_asc_rad)*sin(omega_rad)+cos(right_asc_rad)*cos(omega_rad)*cos(inclination_rad) -cos(right_asc_rad)*sin(inclination_rad)
        sin(omega_rad)*sin(inclination_rad) cos(omega_rad)*sin(inclination_rad) cos(inclination_rad)];
    
    % Convert satellite position from perifocal to ECI
    satellite_eci_pos = orbit_to_ECI * satellite_perifocal_position;
    
    % Convert satellite position from ECI to ECEF for tracking purpose
    ecef_pos = eci2ecef(satellite_eci_pos, t);
    
    sat_ecef(1,position_number) = ecef_pos(1,:);
    sat_ecef(2,position_number) = ecef_pos(2,:);
    sat_ecef(3,position_number) = ecef_pos(3,:);
    
    % Update position number
    position_number = position_number + 1;
    
end

%% Plot ground trace

satellite_llhgd = ecef2llhgd(sat_ecef);

figure
map = imread('Map.jpg');
image([-180 180], [90 -90], map);
axis xy
hold on
lat = rad2deg(satellite_llhgd(1,:));
lon = rad2deg(satellite_llhgd(2,:));

scatter(lon,lat,0.1,'r');

end







