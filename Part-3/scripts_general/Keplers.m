%% function Kepler equations
%
% Author Kuan Chun Hwang

function [pos_eci, pos_llhgc, pos_ecef, theta, semi_rect] = Keplers(n_sat,a,e,i_rad,cap_omega_rad,omega_rad,M_rad,Epoc_time,times,t_vernal)

% Initialize global vairiables
global mu_earth;

n_times = length(times);
n_pos = 3;
sat_numb = length(n_sat);
pos_eci = zeros(n_times, sat_numb, n_pos);
pos_llhgc = zeros(n_times, sat_numb, n_pos);
pos_ecef = zeros(n_times, sat_numb, n_pos);

for n = 1:sat_numb
    
    % Extract Data
    a_n = a(n);
    e_n = e(n);
    i_radn = i_rad(n);
    cap_omega_radn = cap_omega_rad(n);
    M_radn= M_rad(n);
    omega_radn = omega_rad(n);
    epo_time = Epoc_time(n);
    
    % Initialize first guess of eccentric anomaly which is t0 = t so mean
    % anomaly is the same as eccentric anomaly
    e_anomaly_old_rad = M_radn;
    M0 = M_radn;
    
    % Set tolerance which is the accuracy of the iterations and assume
    % difference is great at the start
    toler = 0.00000001;
    
    % Initialize position number
    position_number = 1;
    
    % Orbit plotting loop
    for t=times
        
        % reset delta value
        delta = 10000000000000;
        
        % Set mean anomay at epoch time in that time step
        M = M0 + (sqrt(mu_earth./(a_n.^3))).*(t-epo_time);
        
        % Start Iteration
        while delta > toler
            
            % Solve Keplers equation for E
            delta = e_anomaly_old_rad - e_n .* sin(e_anomaly_old_rad) - M;
            
            % Next iteration value
            e_anomaly_new_rad = e_anomaly_old_rad - delta ./ (1 - e_n .* cos(e_anomaly_old_rad));
            e_anomaly_old_rad = e_anomaly_new_rad;
            
        end
        
        % Final eccentric anomaly value
        e_anomaly_ans = e_anomaly_old_rad;
        
        % Find True Anomaly
        true_anomaly_rad = 2 .* atan2( (sqrt((1+e_n)./(1-e_n)) * tan (e_anomaly_ans./2)),1);
        
        % Find semilatus rectum (p)
        p = a_n .* (1 - e_n.^2);
        
        if t == times(1)
            theta = true_anomaly_rad;
            semi_rect = p;
        end
        
        % Calculate Distance to satellite
        distance = p ./ (1 + e_n .* cos (true_anomaly_rad));
        
        % Coordinates in orbital frame
        x_orbital = distance * cos(true_anomaly_rad);
        y_orbital = distance * sin(true_anomaly_rad);
        z_orbital = 0;
        
        satellite_perifocal_position = [x_orbital; y_orbital; z_orbital];
        
        % Initialize Transfer matrix from orbital coordinates to ECI frame
        orbit_to_ECI = orbit2ECI(cap_omega_radn,omega_radn,i_radn);
        
        % Convert satellite position from perifocal to ECI
        satellite_eci_pos = orbit_to_ECI * satellite_perifocal_position;
        
        % Extract ECI data
        x_eci(position_number) = satellite_eci_pos(1,:);
        y_eci(position_number) = satellite_eci_pos(2,:);
        z_eci(position_number) = satellite_eci_pos(3,:);
        
        % Convert satellite position from ECI to ECEF for tracking purpose
        ecef_pos = eci2ecef(satellite_eci_pos, (t-t_vernal));
        
        % Extrace ecef data
        x_ecef(position_number) = ecef_pos(1,:);
        y_ecef(position_number) = ecef_pos(2,:);
        z_ecef(position_number) = ecef_pos(3,:);
        
        % Convert from ecef to llhgc
        llhgc_pos = ecef2llhgc(ecef_pos);
        
        % Extract llhgc data
        lat(position_number) = rad2deg(llhgc_pos(1,:));
        lon(position_number) = rad2deg(llhgc_pos(2,:));
        alt(position_number) = llhgc_pos(3,:);
        
        % Update position number
        position_number = position_number + 1;
    end
    pos_ecef(:, n, :) = reshape([x_ecef; y_ecef; z_ecef]', n_times, 1, n_pos);
    
    pos_eci(:, n, :) = reshape([ x_eci; y_eci; z_eci]', n_times, 1, n_pos);
    
    pos_llhgc(:, n, :) = reshape([ lat; lon; alt]', n_times, 1, n_pos);
    
end

end