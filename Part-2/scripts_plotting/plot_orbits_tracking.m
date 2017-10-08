%% function 3D orbit plot and ground trace with estimated parameters
%
% Author Kuan Chun Hwang

%% Plot 3D orbit
function plot_orbits_tracking(pos_eci_i,pos_ecef_i,pos_llhgc_i,times,t_stamp,pos_eci_f,pos_ecef_f,pos_llhgc_f,gs_llh)
global w_earth;
% Convert from radians to degree
gs_lat_deg = rad2deg(gs_llh(1,:));
gs_lon_deg = rad2deg(gs_llh(2,:));
%% User input
prompt = {'ECI or ECEF'};
dlg_title = 'Coordinate Frame';
answer = inputdlg(prompt,dlg_title);

%% 3D ECI plot

n_times = length(times);

screen1 = [0.0 0.0 0.5 1.0];
screen2 = [0.0 0.0 0.5 1.0];
screen3 = [0.5 0.0 0.5 0.5];
tr = strcmp('ECI',answer);
if tr == 1
    % Create figure and load topographical Earth map
    fig.globe = figure;
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
    title('3D ECI plot');
    xlabel('x(m)')
    ylabel('y(m)')
    zlabel('z(m)')
    hold on
    sh=surface(x,y,z,props);
    
    set(fig.globe, 'Units', 'normalized', 'Position', screen1);
    plt.orbits = plot3(pos_eci_i(:, :, 1), pos_eci_i(:, :, 2), pos_eci_i(:, :, 3),'r');
    hold on;
    plt.sats = scatter3(NaN, NaN, NaN,'r', 'filled', ...
        'XDataSource', 'pos_eci_i(i, :, 1)', ...
        'YDataSource', 'pos_eci_i(i, :, 2)', ...
        'ZDataSource', 'pos_eci_i(i, :, 3)');
    hold on;
    plt.orbits_est = plot3(pos_eci_f(:, :, 1), pos_eci_f(:, :, 2), pos_eci_f(:, :, 3),'b');
    hold on;
    plt.sats_est = scatter3(NaN, NaN, NaN,'b', 'filled', ...
        'XDataSource', 'pos_eci_f(i, :, 1)', ...
        'YDataSource', 'pos_eci_f(i, :, 2)', ...
        'ZDataSource', 'pos_eci_f(i, :, 3)');
    
    %% 3D ECEF plot
else
    
    % Create figure and load topographical Earth map
    fig.globe = figure;
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
    title('3D ECEF plot');
    xlabel('x(m)')
    ylabel('y(m)')
    zlabel('z(m)')
    hold on
    sph=surface(x,y,z,props);
    
    set(fig.globe, 'Units', 'normalized', 'Position', screen2);
    plt.orbitsECEF = plot3(pos_ecef_i(:, :, 1), pos_ecef_i(:, :, 2), pos_ecef_i(:, :, 3),'r');
    hold on;
    plt.satsECEF = scatter3(NaN, NaN, NaN,'r', 'filled', ...
        'XDataSource', 'pos_ecef_i(i, :, 1)', ...
        'YDataSource', 'pos_ecef_i(i, :, 2)', ...
        'ZDataSource', 'pos_ecef_i(i, :, 3)');
    hold on;
    plt.orbitsECEF_est = plot3(pos_ecef_f(:, :, 1), pos_ecef_f(:, :, 2), pos_ecef_f(:, :, 3),'b');
    hold on;
    plt.satsECEF_est = scatter3(NaN, NaN, NaN,'b', 'filled', ...
        'XDataSource', 'pos_ecef_f(i, :, 1)', ...
        'YDataSource', 'pos_ecef_f(i, :, 2)', ...
        'ZDataSource', 'pos_ecef_f(i, :, 3)');
end
%% 2D ground trace
fig.map = figure;
map = imread('Map.jpg');
image([-180 180], [90 -90], map);
axis xy
hold on
% Plot Ground Station Locations
scatter(gs_lon_deg,gs_lat_deg,'m','LineWIdth',2)
hold on
set(fig.map, 'Units', 'normalized', 'Position', screen3);
plt.ground = plot(pos_llhgc_i(:, :, 2), pos_llhgc_i(:, :, 1),'r');
plt.ground_trace = scatter(NaN, NaN,'r', 'filled', ...
    'XDataSource', 'pos_llhgc_i(i, :, 2)', ...
    'YDataSource', 'pos_llhgc_i(i, :, 1)');
hold on
plt.ground_est = plot(pos_llhgc_f(:, :, 2), pos_llhgc_f(:, :, 1),'b');
plt.ground_trace_est = scatter(NaN, NaN,'b', 'filled', ...
    'XDataSource', 'pos_llhgc_f(i, :, 2)', ...
    'YDataSource', 'pos_llhgc_f(i, :, 1)');
axis equal;
xlim([-180 180]);
ylim([-90 90]);
title('2D Ground Trace plot');
xlabel('longitude (deg)')
ylabel('latitude (deg)')


%% Plot orbits
speedup = 5;
for i = 1:speedup:n_times
    
    if tr == 1
        % Rotate Earth
        rotate(sh,[0,0,1],(speedup*t_stamp*w_earth*180/pi));
        refreshdata(plt.sats, 'caller');
        refreshdata(plt.sats_est, 'caller');
    else
        refreshdata(plt.satsECEF, 'caller');
        refreshdata(plt.satsECEF_est, 'caller');
    end
    refreshdata(plt.ground_trace, 'caller');
    refreshdata(plt.ground_trace_est, 'caller');
    drawnow();
end