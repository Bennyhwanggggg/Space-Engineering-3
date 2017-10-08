%% Part 3
%
% Author Kuan Chun Hwang

clear;
clc;
close all
addpath('./scripts_general','./module_conversion');
warning('off','all');

% Initialize global constants
constants();
global mu_earth;
global v_initial;
global v_final;
global a_i;
global a_f;
global i_f_rad;
global RAAN_f_rad;
global rDot_final;
global r_final;

%% User input
[FDS, CDS, SOD, BFGS, i_desired, RAAN_input, vg, bg, amer, bmer] = userInput();
tic
%% Initialize satellite park orbit
a_i = 6655937;      % Semi-major axis (m)
e_i = 0;            % Eccentricity
i_i_rad = deg2rad(-28.5);    % Inclination (deg)
RAAN_i_rad = 0;     % Right Acension of Ascending Node (rad)
omega_i_rad = 0;    % Argument of Perigee (rad)
M0_i_rad = 0;       % Mean Anomaly (rad)
t0_i = 0;           % Epoch (s)

Initial_paramters = [a_i, e_i, i_i_rad, RAAN_i_rad, omega_i_rad, M0_i_rad,t0_i];

% Calculate initial orbit velocity in satellite parking orbit
v_initial = sqrt(mu_earth/a_i);

%% Desired orbit paramter
% The semi-major axis of the desired final orbit gives the radius required
% for a GEO orbit. Eccentricity is set to 0 to ensure circular orbit while
% inclination and RAAn are controlled variable. Argument of perigee, Mean
% anomaly and epoch are free parameters.
a_f = 42164000;         % Semi-major axis (m)
e_f = 0;                % Eccentricity
i_f_rad = i_desired;    % Inclination (rad)
RAAN_f_rad = RAAN_input;% Right Acension of Ascending Node (rad)
omega_f_rad = 0;        % Argument of Perigee (rad)
M0_f_rad = 0;           % Mean Anomaly (rad)
t0_f = 0;               % Epoch (s)

fprintf('Radius required for a GEO orbit = %g\n',a_f);

Final_paramters = [a_f, e_f, i_f_rad, RAAN_f_rad, omega_f_rad, M0_f_rad,t0_f];

% Calculate final orbit velocity,the velocity required for GEO orbit
v_final = sqrt(mu_earth/a_f);
fprintf('Velocity required for a GEO orbit = %g\n',v_final);

%% Set times
% Calculate the period of each orbit
period_initial = floor(calc_period(a_i));
period_final = floor(calc_period(a_f));
t_stamp = 300;
times_i = 0:t_stamp:period_initial;
times_f = 0:t_stamp:period_final;

%% Calculate Initial and Final orbit ECI positions using known orbital paramters
% Using Keplers equations we find the ECI positions
% Use Kepler equation to find positions in orbit
n_sat = 1;
t_vernal = 0;
[pos_eci_initial, pos_llhgc_initial, pos_ecef_initial, true_anomaly_rad_initial, p_initial] = Keplers(n_sat,a_i,e_i,i_i_rad,RAAN_i_rad,omega_i_rad,M0_i_rad,t0_i,times_i,t_vernal);
[pos_eci_final, pos_llhgc_final, pos_ecef_final, true_anomaly_rad_final, p_final] = Keplers(n_sat,a_f,e_f,-i_f_rad,RAAN_f_rad,omega_f_rad,M0_f_rad,t0_f,times_f,t_vernal);

%% Find position and velocity vector of the initial and final orbit
% Find initial position vector and velocity vector that can be used to
% initiate dynamic model
[r_initial, rDot_initial] = find_r_rdot_ECI(pos_eci_initial,e_i,i_i_rad,omega_i_rad,RAAN_i_rad,true_anomaly_rad_initial,p_initial);
% Find the final orbit vector and velocity vector so they can be used as
% contraints
equinoctialStates = Classical2Equinoctial(a_f, e_f, i_f_rad, omega_f_rad, RAAN_f_rad, pi);
[r_final, rDot_final] = pos_velEquin2ECI(equinoctialStates);
%% Estimate deltaV required for transfer orbit
% Using Hohmann transfer equation (This only applies to 2D so it is a rough
% estimate)
dV1est = sqrt(mu_earth/a_i)*(sqrt((2*a_f)/(a_i+a_f))-1);
dV2est = sqrt(mu_earth/a_f)*(1-sqrt((2*a_i)/(a_i+a_f)));

%% Initialize optimisation parameters
% Use very good estimate or bad esimate depending on user input
if vg == 1
    theta1 = 0;         % Guess first burn elevation angle (rad)
    psi1 = 0;           % Guess first burn azimuth angle (rad)
    theta2 = 0;         % Guess second burn elevation angle (rad)
    psi2 = 0;           % Guess second burn azimuth angle (rad)
    dV1 = dV1est;       % Guess for initial magnitude of first burn velocity change
    dV2 = dV2est;       % Guess for initial magnitude of second burn velocity
    r0  = r_initial;    % Satellite parking orbit r constraint
    v0  = rDot_initial; % Satellite parking orbit v constraint
    dE1 = pi;           % Guess change in eccentric anomaly of initial orbit
    dE2 = pi;           % Guess change in eccentric anomaly of final orbit
elseif bg == 1
    theta1 = 40;        % Guess first burn elevation angle (rad)
    psi1 = 30;          % Guess first burn azimuth angle (rad)
    theta2 = 20;        % Guess second burn elevation angle (rad)
    psi2 = 10;          % Guess second burn azimuth angle (rad)
    dV1 = dV1est*5;     % Guess for initial magnitude of first burn velocity change
    dV2 = dV2est*5;     % Guess for initial magnitude of second burn velocity
    r0  = r_initial;    % Satellite parking orbit r constraint
    v0  = rDot_initial; % Satellite parking orbit v constraint
    dE1 = pi;           % Guess change in eccentric anomaly of initial orbit
    dE2 = pi;           % Guess change in eccentric anomaly of final orbit
end

% Introduce scaling factor
scaleFactor = [1/pi; 1/1000; 1/pi; 1/pi; 1/pi; 1/1000; 1/pi; 1/pi];
x_parameter = scaleFactor.*[dE1, dV1, theta1, psi1, dE2, dV2, theta2, psi2]';

%% Sequential Qudartic Programing (SQP)
% To minmise fuel, the total deltaV required for orbit transfer will need
% to minimised. Therefore, we need to use Sequential Quadratic programing
% to find the optimised parameters that minmise delta V
[x_parameter,L_data,c_data,f_data,alpha_data,error_data,x_data,count] = SQP(x_parameter,r0,v0,scaleFactor,FDS, CDS, SOD, BFGS, amer, bmer);

%% Print final results
% Scale back final result
x_result = x_parameter./scaleFactor;

% Extract results into relevant variables
Optimal_dE1 = rad2deg(x_result(1));
Optimal_dV1 = x_result(2);
Optimal_theta1 = rad2deg(x_result(3));
Optimal_psi1 = rad2deg(x_result(4));
Optimal_dE2 = rad2deg(x_result(5));
Optimal_dV2 = x_result(6);
Optimal_theta2 = rad2deg(x_result(7));
Optimal_psi2 = rad2deg(x_result(8));
Total_dV = Optimal_dV1+Optimal_dV2;
% Print results to 3d.p
fprintf('delta E1      =   %0.3f deg\n',Optimal_dE1);
fprintf('delta V1      =   %0.3f m/s\n',Optimal_dV1);
fprintf('theta1        =   %0.3f deg\n',Optimal_theta1);
fprintf('psi1          =   %0.3f deg\n',Optimal_psi1);
fprintf('delta E2      =   %0.3f deg\n',Optimal_dE2);
fprintf('delta V2      =   %0.3f m/s\n',Optimal_dV2);
fprintf('theta2        =   %0.3f deg\n',Optimal_theta2);
fprintf('psi2          =   %0.3f deg\n',Optimal_psi2);
fprintf('Total Delta V =   %0.3f m/s\n',Total_dV);

%% Plot the cruise orbit and target orbit
plotOrbitTransfer(pos_eci_initial,pos_eci_final,Optimal_dE1,r0,v0,Optimal_dE2,Optimal_theta1,Optimal_psi1,Optimal_dV1,Optimal_dV2,Optimal_theta2,Optimal_psi2);

%% Produce relevant plots
% Relevant plots are produced to show SQP optimisation process
figure
plot(count,alpha_data);
title('Step size change during iteration process')
xlabel('iteration number')
ylabel('Step size')
figure
plot(count,error_data);
title('Error change during iteration process')
xlabel('iteration number')
ylabel('Error')
figure
plot(count,f_data);
title('Objective iteration process')
xlabel('iteration number')
ylabel('Objective')
figure
plot(count,L_data);
title('Lagrangian changes during iteration process')
xlabel('iteration number')
ylabel('Lagrangian size')
figure
plot(count,c_data(1,:),count,c_data(2,:),count,c_data(3,:),count,c_data(4,:),count,c_data(5,:));
title('Constraint changes during iteration process')
xlabel('iteration number')
ylabel('Constraint size')
legend('boundary condition 1','boundary condition 2','boundary condition 3','boundary condition 4','boundary condition 5')
figure
plot(count,x_data(1,:),count,x_data(2,:),count,x_data(3,:),count,x_data(4,:),count,x_data(5,:),count,x_data(6,:),count,x_data(7,:),count,x_data(8,:));
title('Optimisation parameter changes during iteration process')
xlabel('iteration number')
ylabel('Scaled optimisation parameter values')
legend('delta E1','delta V1','theta 1','phi 1','delta E2','delta V2','theta 2','phi 2')
toc