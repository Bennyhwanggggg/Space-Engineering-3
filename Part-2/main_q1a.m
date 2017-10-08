%% 3D orbit and ground trace
%
% Author: Kuan Chun Hwang

%% Close everything
clear;
clc;
close all;

% Add file paths
addpath('./scripts_general', './scripts_earth_motion', ...
    './module_conversion', './scripts_orbit_calculations','./sat_data','./scripts_plotting');

% Obtain global constants
constants();
global mu_earth;

%% 3D orbit and ground trace

% Import data
SatData = importdata('GPSsat_ephem.txt');

n_sat = SatData(:,1);           % Satellite number
a = SatData(:,2);               % Semi major axis (m)
e = SatData(:,3);               % Eccentricity
i_deg = SatData(:,4);           % Inclination (deg)
cap_omega_deg = SatData(:,5);   % Right ascension of ascending node (deg)
omega_deg = SatData(:,6);       % Argument of Perigee (deg)
M_deg = SatData(:,7);           % Mean anomaly at Epoch (deg)
Epoc_time = SatData(:,8);       % Epoch Time (s)

% Set time variables
t_step = 10;                   % (s)
times = 1:t_step:12*60*60;
t_vernal = 0;

% Find position in ECI, llhgc and ECEF
[pos_eci, pos_llhgc, pos_ecef] = keplerian(n_sat,a,e,i_deg, cap_omega_deg, omega_deg, M_deg,Epoc_time,times,t_vernal);

% Plot orbits
plot_orbits(pos_eci,pos_ecef,pos_llhgc,times,t_step)

