%% Part1.3
%
%  Author: Kuan Chun Hwang

% Close everything
clear;
clc;
close all;
% Add file paths
addpath('./scripts_general', './scripts_earth_motion', ...
    './module_conversion', './scripts_orbit_calculations');

% Obtain global constants
constants();

%% Choose a ground station and produce range, azimuth and elevation plot for that location

% Satellite's orbit is assumed to have no perturbation
% Initialize orbital parameter of SARAL (the chosen LEO satellite)
inclination_LEO = 98.5412;  % Degrees
right_asc_LEO = 251.8101;   % Right ascension of the ascending node (degrees)
eccentricity_LEO = 0.0000401;
omega_LEO = 050.0426;       % Argument of Perigee (degrees)
mean_anomaly_LEO = 310.0793;% Degrees
mean_motion_LEO = 14.32253629*2*pi/(24*60*60); %rev/day

orig_para = [inclination_LEO;eccentricity_LEO;right_asc_LEO;omega_LEO;mean_anomaly_LEO;mean_motion_LEO];

% The ground station chosen is Sydney
% % Sydney Latitude: 33.8681° S, 151.19° E and altitude is assumed to be
% sea level
gs_lat = deg2rad(-33.8681);   %rad
gs_lon = deg2rad(151.19);   %rad
gs_alt = 0;

gs_llh = [gs_lat; gs_lon; gs_alt];

% Time stamp per step (s)
t_stamp = 1;
sat_eci_pos = kepl_orbit(inclination_LEO, right_asc_LEO, eccentricity_LEO, omega_LEO, mean_anomaly_LEO, mean_motion_LEO,t_stamp);

% Plot Range, az, el
[sat_range, sat_az, sat_el] = plotRazEL(sat_eci_pos,t_stamp,gs_llh);

sat_LGpolar(1,:) = sat_range;
sat_LGpolar(2,:) = sat_az;
sat_LGpolar(3,:) = sat_el;

max_elevation = max(sat_LGpolar(3,:));
fprintf('Maximum Elevation is = %g\n', max_elevation);

%% Part 2, use tracking to find orbit parameter

% Taking observations to use for prediction of orbital parameters
time_gap = 10;
obs_n1 = 64000;
obs_n2 = obs_n1 + time_gap;
obs_n3 = obs_n2 + time_gap;

tracking_calc(obs_n1,obs_n2,obs_n3,t_stamp,sat_LGpolar,gs_llh,orig_para);

    
    
