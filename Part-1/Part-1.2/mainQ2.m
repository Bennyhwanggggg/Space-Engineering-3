%% Part 1.2
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
global mu_earth;


%% Plot perturbated orbital path of LEO satellite - SARAL

% Initialize orbital parameter of SARAL (the chosen LEO satellite)
inclination = 98.5412;  % Degrees
right_asc = 251.8101;   % Right ascension of the ascending node (degrees)
eccentricity = 0.0000401;
omega = 50.0426;       % Argument of Perigee (degrees)
mean_anomaly = 310.0793;% Degrees
mean_motion = 14.32253629*2*pi/(24*60*60); %rev/day
semi_major_axis = nthroot(mu_earth/(mean_motion^2),3); %(m)

% Convert degress to radians
inclination_rad = deg2rad(inclination);
right_asc_rad = deg2rad(right_asc);
omega_rad = deg2rad(omega);
mean_anomaly_rad = deg2rad(mean_anomaly);

% Find true anomally
% Initialize first guess of eccentric anomaly which is t0 = t so mean
% anomaly is the same as eccentric anomaly
eccentric_anomaly_old_rad = mean_anomaly_rad;

% Set tolerance which is the accuracy of the iterations and assume
% difference is great at the start
toler = 0.0000000000000001;

% Set time step
t_stamp = 50;

% Solve Keplers equation for E
delta = eccentric_anomaly_old_rad - eccentricity * sin(eccentric_anomaly_old_rad) - mean_anomaly_rad;

% true anomaly position
n = 1;

% Start Iteration
for t=1:t_stamp:86400
    
    % Set mean anomay at any given time in that time step
    mean_anomaly_rad = mean_anomaly_rad + (sqrt(mu_earth/(semi_major_axis^3)))*t_stamp;
    
    while abs(delta) > toler
        
        % Next iteration value
        eccentric_anomaly_new_rad = eccentric_anomaly_old_rad - delta / (1 - eccentricity * cos(eccentric_anomaly_old_rad));
        eccentric_anomaly_old_rad = eccentric_anomaly_new_rad;
        delta = eccentric_anomaly_old_rad - eccentricity * sin(eccentric_anomaly_old_rad) - mean_anomaly_rad;
        
    end
    
    % reset delta value
    delta = 10000000000000;
    
    % Final eccentric anomaly value
    essentric_anomaly_ans = eccentric_anomaly_old_rad;
    
    % Find True Anomaly and store the data
    true_anomaly_rad(n)= 2 * atan2( sqrt((1+eccentricity)/(1-eccentricity)) * tan (essentric_anomaly_ans/2),1);
    
    n = n+1;
end

% Convert from classical orbital parameters to equinoctial elements model
p = semi_major_axis * (1 - (eccentricity^2));
f = eccentricity * cos (omega_rad + right_asc_rad);
g = eccentricity * sin (omega_rad + right_asc_rad);
h = tan(inclination_rad/2) * cos(right_asc_rad);
k = tan(inclination_rad/2) * sin(right_asc_rad);
L = right_asc_rad + omega_rad + true_anomaly_rad(1);

% Put equinoctial parameters into a matrix
x_parameters = [ p; f; g; h; k; L];

% Find satellite ECI position using runge-kutta method
eciPos = runge_kutta(x_parameters,t_stamp,inclination_rad,right_asc_rad,omega_rad,semi_major_axis,eccentricity,true_anomaly_rad);

%% Plot Ground Trace
ground_trace_ECI(eciPos,t_stamp);
title('ground trace of satellite')

%% Plot 3D orbit
% Create figure and load topographical Earth map
figure
load('topo.mat','topo');

% Create a sphere, make it earth sized (in meters)
[x,y,z] = sphere(50);
x = -x.*6378000;
y = -y.*6378000;
z = z.*6378000;

props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;

% Plot Earth
axes('dataaspectratio',[1 1 1],'visible','on')
hold on
sh=surface(x,y,z,props);

for position_number=1:(86400/t_stamp)
    
    % Rotate Earth
    rotate(sh,[0,0,1],(t_stamp*w_earth*180/pi));
    
    % Plot orbit
    plot3(eciPos(1,position_number),eciPos(2,position_number),eciPos(3,position_number));
    
    % Plot satellite location and remove old position
    satellite = plot3(eciPos(1,position_number),eciPos(2,position_number),eciPos(3,position_number),'ro','MarkerFaceColor','r');
    pause(0.000000000000001);
    delete(satellite);
    
    % Label graph
    title('3D orbit of satellite')
    xlabel('x (m)') % x-axis label
    ylabel('y (m)') % y-axis label
    zlabel('z (m)') % z-axis label
end
