%% UAV
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
global r_earth;

%% Ground Station Data
% Ground Station Data
lat = deg2rad(-34.76);  %rad
lon = deg2rad(150.03);  %rad
alt = 680;              %(m)
gs_llh = [lat;lon;alt];

% Convert from llh to LGCV
gs_ecef = llhgc2ecef(gs_llh);

%% Satellite positions
% Import GPSsat_ephem data
ephemData = importdata('GPSsat_ephem.txt');

% Extract data
n_sat = ephemData(:,1);           % Satellite number
a = ephemData(:,2);               % Semi major axis (m)
e = ephemData(:,3);               % Eccentricity
i_deg = ephemData(:,4);           % Inclination (deg)
cap_omega_deg = ephemData(:,5);   % Right ascension of ascending node (deg)
omega_deg = ephemData(:,6);       % Argument of Perigee (deg)
M_deg = ephemData(:,7);           % Mean anomaly at Epoch (deg)
Epoc_time = ephemData(:,8);       % Epoch Time (s)

% Import pseudorange data
pseudoData = importdata('GPS_pseudorange_F1.txt');

% Extract data
times_pse = pseudoData(:,1);        % Time (s)
n_sat_pse = pseudoData(:,2);       % Satellite number
pseudor_pse = pseudoData(:,3);     % Pseudorange (m)

% Set t_vernal
t_vernal = 7347737.336;

% Find the Ecef position of the observation satellites
[xsv, ysv, zsv] = find_ecef_obs(n_sat, a, e, i_deg, cap_omega_deg, omega_deg, M_deg,Epoc_time, times_pse, n_sat_pse,t_vernal);

% Find visible satellite polar position w.r.t ground station
vis_sat_ecef = [xsv; ysv; zsv];
vis_sat_translate2gs = bsxfun(@minus, vis_sat_ecef(1:3,:), gs_ecef);
vis_sat_lg_cart_wrt_gs = ecef2lg(vis_sat_translate2gs,gs_llh);

% Convert visible satellite position w.r.t ground station from cartesian to
% polar
vis_sat_lg_polar_wrt_gs = cartesian2polar(vis_sat_lg_cart_wrt_gs);

%% UAV Data
% Import UAV true position data
UAV_trueData = importdata('UAVPosition_F1.txt');

% Extract data
UAVtrue_time = UAV_trueData(:,1);
UAV_x = UAV_trueData(:,2);
UAV_y = UAV_trueData(:,3);
UAV_z = UAV_trueData(:,4);

% Sort cartesian NED into matrix
UAVtrue_lg_cart = [UAV_x'; UAV_y'; UAV_z'];

% Convert from cartesian to polar
UAV_lg_polar = cartesian2polar(UAVtrue_lg_cart);

%% Find estimated position of UAVF1

% Find the times in psudorange data that are in the true data and extract
% those points and find number of visible satellites and time of observations
[rep_t,vis_sat_n,t_obs] = find_t_repeat(times_pse);
% Fix error
t_used = rep_t;
t_used(1) = [];

% Find estimated UAV position, H matrix and the time sets
[x_final,Hs,time_set] = pesudorange(times_pse,pseudor_pse,xsv,ysv,zsv);
% Fix Error
x_final(:,1) = [];
% Extract Clock Bias
ClockBias = x_final(4,:);

%% Find estimated trajectory

% Convert estimated points into LGCV
est_ecef = x_final(1:3,:);
est_ecef_cart = bsxfun(@minus, est_ecef(1:3,:), gs_ecef);
est_lg_cart_wrt_gs = ecef2lg(est_ecef_cart,gs_llh);

% Convert from cartesian coordiantes to polar
est_lg_polar = cartesian2polar(est_lg_cart_wrt_gs);

%% DOPs
for i = 1:length(Hs)
    DOPs(:,i) = findDOPs(Hs{i});
end

%% Analysis Plots

% Polar plots of UAV's trajectory w.r.t from the ground station animation
plotpolars(UAV_lg_polar,est_lg_polar,vis_sat_lg_polar_wrt_gs,time_set)

% Do analysis plots
analysisPlots(t_used,UAVtrue_time,est_lg_cart_wrt_gs,UAVtrue_lg_cart,est_lg_polar,UAV_lg_polar,t_obs,vis_sat_n,ClockBias)

% Plot DOPs and the indexes of time step where these configuration occurs
[maxDOP,minDOP] = plotDOPs(DOPs,rep_t);

% Best and Worst satellite configuration plot
DOPconfig_plot(maxDOP,minDOP,vis_sat_lg_polar_wrt_gs,time_set);








