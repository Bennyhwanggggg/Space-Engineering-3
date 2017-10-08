%% function 3D orbit plot and ground trace
%
% Author Kuan Chun Hwang

%% Plot 3D orbit
function plot_orbits(pos_eci,pos_ecef,pos_llhgc,times,t_stamp)
global w_earth;

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
    plt.orbits = plot3(pos_eci(:, :, 1), pos_eci(:, :, 2), pos_eci(:, :, 3));
    hold on;
    plt.sats = scatter3(NaN, NaN, NaN, 'filled', ...
        'XDataSource', 'pos_eci(i, :, 1)', ...
        'YDataSource', 'pos_eci(i, :, 2)', ...
        'ZDataSource', 'pos_eci(i, :, 3)');
    
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
    plt.orbitsECEF = plot3(pos_ecef(:, :, 1), pos_ecef(:, :, 2), pos_ecef(:, :, 3));
    hold on;
    plt.satsECEF = scatter3(NaN, NaN, NaN, 'filled', ...
        'XDataSource', 'pos_ecef(i, :, 1)', ...
        'YDataSource', 'pos_ecef(i, :, 2)', ...
        'ZDataSource', 'pos_ecef(i, :, 3)');
end
%% 2D ground trace
fig.map = figure;
map = imread('Map.jpg');
image([-180 180], [90 -90], map);
axis xy
hold on

set(fig.map, 'Units', 'normalized', 'Position', screen3);
plt.ground = plot(pos_llhgc(:, :, 2), pos_llhgc(:, :, 1));
plt.ground_trace = scatter(NaN, NaN, 'filled', ...
    'XDataSource', 'pos_llhgc(i, :, 2)', ...
    'YDataSource', 'pos_llhgc(i, :, 1)');
axis equal;
xlim([-180 180]);
ylim([-90 90]);
title('2D Ground Trace plot');
xlabel('longitude (deg)')
ylabel('latitude (deg)')


%% Plot orbits
speedup = 10;
for i = 1:speedup :n_times
    
    if tr == 1
        % Rotate Earth
        rotate(sh,[0,0,1],(speedup*t_stamp*w_earth*180/pi));
        refreshdata(plt.sats, 'caller');
    else
        refreshdata(plt.satsECEF, 'caller');
    end
    refreshdata(plt.ground_trace, 'caller');
    drawnow();
end