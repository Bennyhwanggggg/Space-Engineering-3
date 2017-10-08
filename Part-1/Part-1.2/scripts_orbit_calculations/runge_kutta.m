%% function Runge-Kutta method
%
%  Author: Kuan Chun Hwang

function eciPos = runge_kutta(x_parameters,t_stamp,inclination_rad,right_asc_rad,omega_rad,semi_major_axis,eccentricity,true_anomaly_rad)

global mu_earth;

% Set satellite position number
position_number = 1;

% Start calculating position vector for each time step
for t=1:t_stamp:86400
    
    % Initialize x_k0
    x_k0 = x_parameters;
    
    % Calculate f(x) which is used in the integration part of the equation
    % of motion
    fx_dot = sate_rate_eqn(x_parameters);
    
    % Start Runge-Kutta method
    % Find k1
    k1 = t_stamp * fx_dot;
    
    % Find k2
    k2_parameters = x_k0 + 0.5*k1;
    fx_dotk2 = sate_rate_eqn(k2_parameters);
    k2 = t_stamp * fx_dotk2;
    
    % Find k3
    k3_parameters = x_k0 + 0.5*k2;
    fx_dotk3 = sate_rate_eqn(k3_parameters);
    k3 = t_stamp * fx_dotk3;
    
    % Find k4
    k4_parameters = x_k0 + k3;
    fx_dotk4 = sate_rate_eqn(k4_parameters);
    k4 = t_stamp * fx_dotk4;
    
    % Runge Kutta Method
    x_parameters = x_k0 + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    % Find position vector by converting from Equinoctial elements to ECI position vector
    position_vector = equin2eci(x_parameters);
    
    % Extract Data from position vector    
    eciPos(1,position_number) = position_vector(1);
    eciPos(2,position_number) = position_vector(2);
    eciPos(3,position_number) = position_vector(3);
    
    % Store parameter values for plotting each properties
    p_data(position_number) = x_parameters(1);
    f_data(position_number) = x_parameters(2);
    g_data(position_number) = x_parameters(3);
    h_data(position_number) = x_parameters(4);
    k_data(position_number) = x_parameters(5);
    L_data(position_number) = x_parameters(6);
    
    % Convert Equinoctial elements to Classical
    a_data(position_number) = p_data(position_number) / (1- f_data(position_number)^2 - g_data(position_number)^2);
    e_data(position_number) = sqrt(f_data(position_number)^2 + g_data(position_number)^2);
    i_data(position_number) = atan2(2*sqrt(h_data(position_number)^2+k_data(position_number)^2),1-h_data(position_number)^2-k_data(position_number)^2);
    omega_data(position_number) = atan2(g_data(position_number)*h_data(position_number) - f_data(position_number) * k_data(position_number), f_data(position_number)*h_data(position_number) + g_data(position_number)*k_data(position_number));
    cap_omega_data(position_number) = atan2(k_data(position_number),h_data(position_number));
    v_data(position_number)= L_data(position_number) - atan2(g_data(position_number),f_data(position_number));
    
    % Update satellite position
    position_number = position_number + 1;
end
a_avg = (max(a_data)+ min(a_data))/2;
orbital_period = 2*pi*sqrt((a_avg^3)/mu_earth);
fprintf('Orbital period is = %g\n', orbital_period);

t=1:t_stamp:86400;

figure
subplot(3,2,1);
plot(t,a_data,'b',t,semi_major_axis,'r');
xlabel('time (s)')
ylabel('m')
title('Semi Major Axis')

subplot(3,2,2);
plot(t,e_data,'b',t,eccentricity,'r');
xlabel('time (s)')
ylabel('eccentricity')
title('eccentricity')

subplot(3,2,3);
plot(t,i_data,'b',t,inclination_rad,'r');
xlabel('time (s)')
ylabel('rad')
title('inclination')

subplot(3,2,4);
plot(t,omega_data,'b',t,omega_rad,'r');
xlabel('time (s)')
ylabel('rad')
title('Argument of perigee')

subplot(3,2,5);
plot(t,cap_omega_data,'b',t,right_asc_rad,'r');
xlabel('time (s)')
ylabel('rad')
title('Right Ascension of Ascending Node')

subplot(3,2,6);
plot(t,v_data,'b',t,true_anomaly_rad,'r');
xlabel('time (s)')
ylabel('rad')
title('True Anomaly')
end
