%% function 'Ground Trace of satellite using ECI'
%
% Plots the orbit of satellite in ECEF
%
% Input  :   inclination, right_asc, eccentricity, omega,
%            mean_anomaly, mean_motion, semi_major_axis
% Outputs:  3D plot of orbit in ECEF position
%
% Kuan Chun Hwang

%% This orbital calculation function will plot 3D orbit plot and ground tracking plot
function ground_trace_ECI(eciPos,t_stamp)

n = 1;

for t = 1:t_stamp:86400;
    
    % Convert satellite position from ECI to ECEF for tracking purpose
    ecef_pos = eci2ecef(eciPos, t);
    
    sat_ecef(1,n) = ecef_pos(1,n);
    sat_ecef(2,n) = ecef_pos(2,n);
    sat_ecef(3,n) = ecef_pos(3,n);
    
    % Update position number
    n = n + 1;
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







