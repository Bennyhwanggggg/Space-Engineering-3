%% function "Plot range, aimuth and elevation"
%
%  Author: Kuan- Chun Hwang

function [sat_range, sat_az, sat_el] = plotRazEL(sat_eci_pos,t_stamp,gs_llh)

% Find satellite ecef
gs_ecef = llhgd2ecef(gs_llh);

x_obs = 0;
n = 1;

for t = 0:t_stamp:86400
    
    % Convert ground station coordinates from ECEF to ECI
    gs_eci = ecef2eci(gs_ecef,n*t_stamp);
    
    % Caculate satellites relative position to ground station in ECI
    sat_eci_rel = sat_eci_pos(:,n) - gs_eci;
    
    % Convert satellite position relative to ground station from eci to ecef
    sat_ecef = eci2ecef(sat_eci_rel,(n-1)*t_stamp);
    
    % Convert satellite relative position to GS from ECEF to LGDV
    sat_lg_cart = ecef2lg(sat_ecef,gs_llh);
    
    % Convert from cartesian to polar
    sat_lg_polar = cartesian2polar(sat_lg_cart);
    
    % Store data
    sat_range(n) = sat_lg_polar(1,:);
    sat_az(n) = sat_lg_polar(2,:);
    sat_el(n) = sat_lg_polar(3,:);
    
    % Find satellite ECEF position
    sat_ECEF_notrel = eci2ecef(sat_eci_pos(:,n),(n-1)*t_stamp);    
    
    if sat_el(n) >= 0 && sat_el(n) <= pi
        sat_el_obs(n) = sat_el(n);
        sat_el_noObs(n) = NaN;
        satellite_llhgd_obs(:,n) = ecef2llhgd(sat_ECEF_notrel);
        x_obs = x_obs + 1;
    else
        sat_el_obs(n) = NaN;
        sat_el_noObs(n) = sat_el(n);
        satellite_llhgd(:,n) = ecef2llhgd(sat_ECEF_notrel);
    end
    
    
    % Update position number
    n = n +1;
end
fprintf('Percentage of time in view = %g\n', 100*(x_obs*t_stamp)/86400);

% Plot relevant data
t =  0:t_stamp:86400;
figure(1)
subplot(3,1,1);
plot(t,sat_range);
% Label graph
title('Range vs Time')
ylabel('Range (m)') % x-axis label
xlabel('Time (s)') % y-axis label

subplot(3,1,2);
plot(t,sat_az);
% Label graph
title('azimuth vs Time')
ylabel('azimuth(rad)') % x-axis label
xlabel('Time (s)') % y-axis label

subplot(3,1,3);
plot(t,sat_el_obs,'b',t,sat_el_noObs,'r');
% Label graph
title('Elevation vs Time')
ylabel('Elevation(rad)') % x-axis label
xlabel('Time (s)') % y-axis label

%% Plot ground trace

figure
map = imread('Map.jpg');
image([-180 180], [90 -90], map);
axis xy
hold on
lat_obs = rad2deg(satellite_llhgd_obs(1,:));
lon_obs = rad2deg(satellite_llhgd_obs(2,:));
scatter(lon_obs,lat_obs,1,'r');
hold on
lat = rad2deg(satellite_llhgd(1,:));
lon = rad2deg(satellite_llhgd(2,:));
scatter(lon,lat,1,'y');
title('ground trace');
xlabel('longitude');
ylabel('latitude');
end