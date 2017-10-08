%% Analysis Plots - 3D UAV pos
%                   2D UAV pos
%                   Altitude vs time,
%                   Range Azimuth Elevation vs Time
%                   Polar UAV pos w.r.t Ground station
%                   Number of visible satellites vs Time
%
% Author Kuan Chun Hwang

function analysisPlots(t_used,true_time,est_lg_cart_wrt_gs,true_lg_cart,est_lg_polar,true_lg_polar,t_obs,vis_sat_n,ClockBias)

UAVF1_x = true_lg_cart(1,:);
UAVF1_y = true_lg_cart(2,:);
UAVF1_z = true_lg_cart(3,:);

% Plot 3D position w.r.t Ground Station
figure
plot3(est_lg_cart_wrt_gs(1,:),est_lg_cart_wrt_gs(2,:),est_lg_cart_wrt_gs(3,:),'b');
hold on
plot3(UAVF1_x,UAVF1_y,UAVF1_z,'r');
hold on
scatter3(0,0,0,'fill');

title('3D UAV Position w.r.t Ground Station');
xlabel('x(m)');
ylabel('y(m)');
zlabel('z(m)');
legend('Estimated UAV position','true UAV position','Ground Station position');
axis equal

% Plot 2D UAV Position w.r.t Ground Station (xy plane)
figure
plot(est_lg_cart_wrt_gs(1,:),est_lg_cart_wrt_gs(2,:),'b');
hold on
plot(UAVF1_x,UAVF1_y,'r');
hold on
scatter(0,0,'fill');

title('2D UAV Position w.r.t Ground Station (xy plane)');
xlabel('x(m)');
ylabel('y(m)');
legend('Estimated UAV position','true UAV position','Ground Station position');
axis equal

% Plot UAV altitude vs Time w.r.t Ground Station
figure
plot(t_used,est_lg_cart_wrt_gs(3,:),'b')
hold on
plot(true_time,UAVF1_z,'-r')

title('UAV Altitude vs Time w.r.t Ground Station');
xlabel('time (s)');
ylabel('Altitude (m)');
legend('Estimated UAV Altitude','true UAV Altitude');
axis equal

% Range Azimuth Elevation vs Time w.r.t Ground Station
figure
subplot(3,1,1);
plot(true_time, true_lg_polar(1,:),'-r')
hold on
plot(t_used, est_lg_polar(1,:),'b')

title('Range vs time')
xlabel('time')
ylabel('Range (m)')

subplot(3,1,2);
plot(true_time, true_lg_polar(2,:),'-r')
hold on
plot(t_used, est_lg_polar(2,:),'b')

title('azimuth vs time')
xlabel('time')
ylabel('az(rad)')

subplot(3,1,3);
plot(true_time, true_lg_polar(3,:),'-r')
hold on
plot(t_used, est_lg_polar(3,:),'b')

title('elongation vs time')
xlabel('time')
ylabel('el(rad)')
legend('True','Estimated')

figure
plot(t_obs,vis_sat_n)
title('Number of visible satellites vs time')
xlabel('time')
ylabel('Number of satellites')
axis([min(t_obs) max(t_obs) 0 max(vis_sat_n)+2]);

figure
plot(t_used,ClockBias)
title('UAV Clock Bias vs Time w.r.t Ground station')
xlabel('time')
ylabel('Clock Bias')
end