%% function plot tracking measurements
%
% Author Kuan Chun Hwang

function plot_tracking_measurements(pos_eci_i,pos_eci_f,times_f)

% Find average satellite position error
x_true = pos_eci_i(:,1,1);
y_true = pos_eci_i(:,1,2);
z_true = pos_eci_i(:,1,3);
x_est = pos_eci_f(:,1,1);
y_est = pos_eci_f(:,1,2);
z_est = pos_eci_f(:,1,3);
average_err = norm(sqrt((x_true - x_est).^2 + (y_true - y_est).^2 + (z_true - z_est).^2))/length(times_f);
fprintf('Average satellite position error = %g\n', average_err);

figure
subplot(3,1,1)
plot(times_f,x_true - x_est)
xlabel('times (s)')
ylabel('x error (m)')
title('x distance error')
subplot(3,1,2)
plot(times_f,y_true - y_est)
xlabel('times (s)')
ylabel('y error (m)')
title('y distance error')
subplot(3,1,3)
plot(times_f,z_true - z_est)
xlabel('times (s)')
ylabel('z error (m)')
title('z distance error')

figure
subplot(3,1,1)
plot(times_f,x_true,'r',times_f,x_est,'b')
xlabel('times (s)')
ylabel('x(m)')
title('x location')
subplot(3,1,2)
plot(times_f,y_true,'r',times_f,y_est,'b')
xlabel('times (s)')
ylabel('y(m)')
title('y location')
subplot(3,1,3)
plot(times_f,z_true,'r',times_f,z_est,'b')
xlabel('times (s)')
ylabel('z(m)')
title('z location')
legend('true','estimated')