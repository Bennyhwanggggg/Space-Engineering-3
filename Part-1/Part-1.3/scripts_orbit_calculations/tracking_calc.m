%% function plot tracking effects

function tracking_calc(obs_n1,obs_n2,obs_n3,t_stamp,sat_LGpolar,gs_llh,orig_para)

% Find error when signal to noise ratio is 200db
noise_test = 200;
% Input white Gaussian noise
sat_LG_noise = awgn(sat_LGpolar,noise_test,'measured');
orb_parameters = Herric_Gibbs(obs_n1,obs_n2,obs_n3,t_stamp,sat_LG_noise,gs_llh);

% Caculate the error
differnce = finderror(orig_para,orb_parameters);

disp_error(differnce);

% Add noise variation
nt = 1;
min_db = -50;
max_db = 400;
for noise = min_db:10:max_db;
    
    % Input white Gaussian noise
    sat_LG_noise = awgn(sat_LGpolar,noise,'measured');
    orb_parameters = Herric_Gibbs(obs_n1,obs_n2,obs_n3,t_stamp,sat_LG_noise,gs_llh);
    
    % Caculate the error
    diff = finderror(orig_para,orb_parameters);
    
    % Extract error data
    a_error(nt) = diff(1,:);
    i_error(nt) = diff(2,:);
    e_error(nt) = diff(3,:);
    cap_omega_error(nt) = diff(4,:);
    omega_error(nt) = diff(5,:);
    M_error(nt) = diff(6,:);
    
    nt = nt + 1;
end
%% Plot paramter error with noise variation
noise = min_db:10:max_db;
figure
subplot(2,3,1);
plot(noise,a_error);
xlabel('db');
ylabel('%');
title('Noise affect on Semi Major Axis tracking');

subplot(2,3,2);
plot(noise,i_error);
xlabel('db');
ylabel('%');
title('Noise affect on inclination tracking');

subplot(2,3,3);
plot(noise,e_error);
xlabel('db');
ylabel('e');
title('Noise affect on eccentricity');

subplot(2,3,4);
plot(noise,cap_omega_error);
xlabel('db');
ylabel('%');
title('Noise affect on right ascension of the ascending node');

subplot(2,3,5);
plot(noise,omega_error);
xlabel('db');
ylabel('%');
title('Noise affect on argument of perigee');

subplot(2,3,6);
plot(noise,M_error);
xlabel('db');
ylabel('%');
title('Noise affect on mean anomaly');