%% AERO4701 Assignment 1 Question 1
%
%  Author: Kuan Chun Hwang
%%
% Close everything
clear;
clc;
close all;

% Add file paths
addpath('./scripts_general', './scripts_earth_motion', ...
        './module_conversion', './scripts_orbit_calculations');
    
% Obtain global constants
constants();

%% Plot Orbital path of LEO satellite - SARAL

% Initialize orbital parameter of SARAL (the chosen LEO satellite)
inclination_LEO = 98.5412;  % Degrees
right_asc_LEO = 251.8101;   % Right ascension of the ascending node (degrees)
eccentricity_LEO = 0.0000401;
omega_LEO = 050.0426;       % Argument of Perigee (degrees)
mean_anomaly_LEO = 310.0793;% Degrees
mean_motion_LEO = 14.32253629*2*pi/(24*60*60); %rev/day

% Time stamp per step (s)
t_stamp = 1;

% Plot ground trace of SARAL's orbit
ground_trace(inclination_LEO, right_asc_LEO, eccentricity_LEO, omega_LEO, mean_anomaly_LEO, mean_motion_LEO,t_stamp);
title('Ground trace of LEO satellite - SARAL')

% Increase time stamp per step to reduce running time of the code
t_stamp = 60;

% Plot SARAL's orbit in 3D
keplerian_orbital_calc(inclination_LEO, right_asc_LEO, eccentricity_LEO, omega_LEO, mean_anomaly_LEO, mean_motion_LEO,t_stamp);
title('3D orbit of LEO satellite - SARAL')
%% Plot Orbital path of the MEO satellite - O3B FM07

% Initialize orbital parameter of O3B FM07 (the chosen MEO satellite)
inclination_MEO = 0.0359;  % Degrees
right_asc_MEO = 353.2540;   % Right ascension of the ascending node (degrees)
eccentricity_MEO = 0.0002445;
omega_MEO = 324.0770 ;       % Argument of Perigee (degrees)
mean_anomaly_MEO = 42.6465;% Degrees
mean_motion_MEO = 05.0011571*2*pi/(24*60*60); %rev/day


% Time stamp per step (s)
t_stamp = 1;

% Plot ground trace of O3B FM07's orbit
ground_trace(inclination_MEO, right_asc_MEO, eccentricity_MEO, omega_MEO, mean_anomaly_MEO, mean_motion_MEO,t_stamp);
title('Ground trace of MEO satellite - O3B FM07')

% Increase time stamp per step to reduce running time of the code
t_stamp = 60;

% Plot O3B FM07's orbit in 3D
keplerian_orbital_calc(inclination_MEO, right_asc_MEO, eccentricity_MEO, omega_MEO, mean_anomaly_MEO, mean_motion_MEO,t_stamp);
title('3D orbit of MEO satellite - O3B FM07')
