%% function 3D orbit plot and ground trace with ground station location and visiual lines
%
% Author Kuan Chun Hwang

%% Plot 3D orbit
function plot_orbits_with_gs(pos_eci,pos_ecef,pos_llhgc,times,t_stamp,gs_llh,visible_time)
global w_earth;

% Convert from radians to degree
gs_lat_deg = rad2deg(gs_llh(1,:));
gs_lon_deg = rad2deg(gs_llh(2,:));

% Convert from llhgc to ecef
gs_ecef = llhgc2ecef(gs_llh);

n_times_gs = length(times);
gs_num = length(gs_ecef);
pos_n = 3;
gs_eci = zeros(n_times_gs,gs_num ,pos_n);

% Convert from ecef to eci
for gs_number = 1:length(gs_ecef)
    for i = 1:length(times)
        grounds_eci(:,i) = ecef2eci(gs_ecef(:,gs_number),times(i));
    end
    gs_eci(:, gs_number, :) = reshape(grounds_eci', n_times_gs, 1, pos_n);
end
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
    hold on;
    plt.gs = scatter3(NaN, NaN, NaN, 'filled', ...
        'XDataSource', 'gs_eci(i, :, 1)', ...
        'YDataSource', 'gs_eci(i, :, 2)', ...
        'ZDataSource', 'gs_eci(i, :, 3)');
    
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
    
    hold on
    scatter3(gs_ecef(1,:),gs_ecef(2,:),gs_ecef(3,:));
    
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
% Plot Ground Station Locations
scatter(gs_lon_deg,gs_lat_deg,'r','LineWIdth',2)

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
speedup = 20;
for i = 1:speedup:n_times
    
    if tr == 1
        % Rotate Earth
        rotate(sh,[0,0,1],(speedup*t_stamp*w_earth*180/pi));
        refreshdata(plt.sats, 'caller');
        refreshdata(plt.gs, 'caller');
        
        % GS 1
        visT1_1 =  visible_time{1,1};
        if any(abs(visT1_1-i)<0.0001)
            visL1_1 = quiver3(pos_eci(i, 1, 1),pos_eci(i, 1, 2),pos_eci(i, 1, 3),gs_eci(i, 1, 1)-pos_eci(i, 1, 1),gs_eci(i, 1, 2)-pos_eci(i, 1, 2),gs_eci(i, 1, 3)-pos_eci(i, 1, 3),'r');
        end
        visT1_2 =  visible_time{1,2};
        if any(abs(visT1_2-i)<0.0001)
            visL1_2 = quiver3(pos_eci(i, 2, 1),pos_eci(i, 2, 2),pos_eci(i, 2, 3),gs_eci(i, 1, 1)-pos_eci(i, 2, 1),gs_eci(i, 1, 2)-pos_eci(i, 2, 2),gs_eci(i, 1, 3)-pos_eci(i, 2, 3),'r');
        end
        visT1_3 =  visible_time{1,3};
        if any(abs(visT1_3-i)<0.0001)
            visL1_3 = quiver3(pos_eci(i, 3, 1),pos_eci(i, 3, 2),pos_eci(i, 3, 3),gs_eci(i, 1, 1)-pos_eci(i, 3, 1),gs_eci(i, 1, 2)-pos_eci(i, 3, 2),gs_eci(i, 1, 3)-pos_eci(i, 3, 3),'r');
        end
        visT1_4 =  visible_time{1,4};
        if any(abs(visT1_4-i)<0.0001)
            visL1_4 = quiver3(pos_eci(i, 4, 1),pos_eci(i, 4, 2),pos_eci(i, 4, 3),gs_eci(i, 1, 1)-pos_eci(i, 4, 1),gs_eci(i, 1, 2)-pos_eci(i, 4, 2),gs_eci(i, 1, 3)-pos_eci(i, 4, 3),'r');
        end
        visT1_5 =  visible_time{1,5};
        if any(abs(visT1_5-i)<0.0001)
            visL1_5 = quiver3(pos_eci(i, 5, 1),pos_eci(i, 5, 2),pos_eci(i, 5, 3),gs_eci(i, 1, 1)-pos_eci(i, 5, 1),gs_eci(i, 1, 2)-pos_eci(i, 5, 2),gs_eci(i, 1, 3)-pos_eci(i, 5, 3),'r');
        end
        visT1_6 =  visible_time{1,6};
        if any(abs(visT1_6-i)<0.0001)
            visL1_6 = quiver3(pos_eci(i, 6, 1),pos_eci(i, 6, 2),pos_eci(i, 6, 3),gs_eci(i, 1, 1)-pos_eci(i, 6, 1),gs_eci(i, 1, 2)-pos_eci(i, 6, 2),gs_eci(i, 1, 3)-pos_eci(i, 6, 3),'r');
        end
        visT1_7 =  visible_time{1,7};
        if any(abs(visT1_7-i)<0.0001)
            visL1_7 = quiver3(pos_eci(i, 7, 1),pos_eci(i, 7, 2),pos_eci(i, 7, 3),gs_eci(i, 1, 1)-pos_eci(i, 7, 1),gs_eci(i, 1, 2)-pos_eci(i, 7, 2),gs_eci(i, 1, 3)-pos_eci(i, 7, 3),'r');
        end
        visT1_8 =  visible_time{1,8};
        if any(abs(visT1_8-i)<0.0001)
            visL1_8 = quiver3(pos_eci(i, 8, 1),pos_eci(i, 8, 2),pos_eci(i, 8, 3),gs_eci(i, 1, 1)-pos_eci(i, 8, 1),gs_eci(i, 1, 2)-pos_eci(i, 8, 2),gs_eci(i, 1, 3)-pos_eci(i, 8, 3),'r');
        end
        visT1_9 =  visible_time{1,9};
        if any(abs(visT1_9-i)<0.0001)
            visL1_9 = quiver3(pos_eci(i, 9, 1),pos_eci(i, 9, 2),pos_eci(i, 9, 3),gs_eci(i, 1, 1)-pos_eci(i, 9, 1),gs_eci(i, 1, 2)-pos_eci(i, 9, 2),gs_eci(i, 1, 3)-pos_eci(i, 9, 3),'r');
        end
        visT1_10 =  visible_time{1,10};
        if any(abs(visT1_10-i)<0.0001)
            visL1_10 = quiver3(pos_eci(i, 10, 1),pos_eci(i, 10, 2),pos_eci(i, 10, 3),gs_eci(i, 1, 1)-pos_eci(i, 10, 1),gs_eci(i, 1, 2)-pos_eci(i, 10, 2),gs_eci(i, 1, 3)-pos_eci(i, 10, 3),'r');
        end
        visT1_11 =  visible_time{1,11};
        if any(abs(visT1_11-i)<0.0001)
            visL1_11 = quiver3(pos_eci(i, 11, 1),pos_eci(i, 11, 2),pos_eci(i, 11, 3),gs_eci(i, 1, 1)-pos_eci(i, 11, 1),gs_eci(i, 1, 2)-pos_eci(i, 11, 2),gs_eci(i, 1, 3)-pos_eci(i, 11, 3),'r');
        end
        visT1_12 =  visible_time{1,12};
        if any(abs(visT1_12-i)<0.0001)
            visL1_12 = quiver3(pos_eci(i, 12, 1),pos_eci(i, 12, 2),pos_eci(i, 12, 3),gs_eci(i, 1, 1)-pos_eci(i, 12, 1),gs_eci(i, 1, 2)-pos_eci(i, 12, 2),gs_eci(i, 1, 3)-pos_eci(i, 12, 3),'r');
        end
        
        % GS 2
        visT2_1 =  visible_time{2,1};
        if any(abs(visT2_1-i)<0.0001)
            visL2_1 = quiver3(pos_eci(i, 1, 1),pos_eci(i, 1, 2),pos_eci(i, 1, 3),gs_eci(i, 2, 1)-pos_eci(i, 1, 1),gs_eci(i, 2, 2)-pos_eci(i, 1, 2),gs_eci(i, 2, 3)-pos_eci(i, 1, 3),'r');
        end
        visT2_2 =  visible_time{2,2};
        if any(abs(visT2_2-i)<0.0001)
            visL2_2 = quiver3(pos_eci(i, 2, 1),pos_eci(i, 2, 2),pos_eci(i, 2, 3),gs_eci(i, 2, 1)-pos_eci(i, 2, 1),gs_eci(i, 2, 2)-pos_eci(i, 2, 2),gs_eci(i, 2, 3)-pos_eci(i, 2, 3),'r');
        end
        visT2_3 =  visible_time{2,3};
        if any(abs(visT2_3-i)<0.0001)
            visL2_3 = quiver3(pos_eci(i, 3, 1),pos_eci(i, 3, 2),pos_eci(i, 3, 3),gs_eci(i, 2, 1)-pos_eci(i, 3, 1),gs_eci(i, 2, 2)-pos_eci(i, 3, 2),gs_eci(i, 2, 3)-pos_eci(i, 3, 3),'r');
        end
        visT2_4 =  visible_time{2,4};
        if any(abs(visT2_4-i)<0.0001)
            visL2_4 = quiver3(pos_eci(i, 4, 1),pos_eci(i, 4, 2),pos_eci(i, 4, 3),gs_eci(i, 2, 1)-pos_eci(i, 4, 1),gs_eci(i, 2, 2)-pos_eci(i, 4, 2),gs_eci(i, 2, 3)-pos_eci(i, 4, 3),'r');
        end
        visT2_5 =  visible_time{2,5};
        if any(abs(visT2_5-i)<0.0001)
            visL2_5 = quiver3(pos_eci(i, 5, 1),pos_eci(i, 5, 2),pos_eci(i, 5, 3),gs_eci(i, 2, 1)-pos_eci(i, 5, 1),gs_eci(i, 2, 2)-pos_eci(i, 5, 2),gs_eci(i, 2, 3)-pos_eci(i, 5, 3),'r');
        end
        visT2_6 =  visible_time{2,6};
        if any(abs(visT2_6-i)<0.0001)
            visL2_6 = quiver3(pos_eci(i, 6, 1),pos_eci(i, 6, 2),pos_eci(i, 6, 3),gs_eci(i, 2, 1)-pos_eci(i, 6, 1),gs_eci(i, 2, 2)-pos_eci(i, 6, 2),gs_eci(i, 2, 3)-pos_eci(i, 6, 3),'r');
        end
        visT2_7 =  visible_time{2,7};
        if any(abs(visT2_7-i)<0.0001)
            visL2_7 = quiver3(pos_eci(i, 7, 1),pos_eci(i, 7, 2),pos_eci(i, 7, 3),gs_eci(i, 2, 1)-pos_eci(i, 7, 1),gs_eci(i, 2, 2)-pos_eci(i, 7, 2),gs_eci(i, 2, 3)-pos_eci(i, 7, 3),'r');
        end
        visT2_8 =  visible_time{2,8};
        if any(abs(visT2_8-i)<0.0001)
            visL2_8 = quiver3(pos_eci(i, 8, 1),pos_eci(i, 8, 2),pos_eci(i, 8, 3),gs_eci(i, 2, 1)-pos_eci(i, 8, 1),gs_eci(i, 2, 2)-pos_eci(i, 8, 2),gs_eci(i, 2, 3)-pos_eci(i, 8, 3),'r');
        end
        visT2_9 =  visible_time{2,9};
        if any(abs(visT2_9-i)<0.0001)
            visL2_9 = quiver3(pos_eci(i, 9, 1),pos_eci(i, 9, 2),pos_eci(i, 9, 3),gs_eci(i, 2, 1)-pos_eci(i, 9, 1),gs_eci(i, 2, 2)-pos_eci(i, 9, 2),gs_eci(i, 2, 3)-pos_eci(i, 9, 3),'r');
        end
        visT2_10 =  visible_time{2,10};
        if any(abs(visT2_10-i)<0.0001)
            visL2_10 = quiver3(pos_eci(i, 10, 1),pos_eci(i, 10, 2),pos_eci(i, 10, 3),gs_eci(i, 2, 1)-pos_eci(i, 10, 1),gs_eci(i, 2, 2)-pos_eci(i, 10, 2),gs_eci(i, 2, 3)-pos_eci(i, 10, 3),'r');
        end
        visT2_11 =  visible_time{2,11};
        if any(abs(visT2_11-i)<0.0001)
            visL2_11 = quiver3(pos_eci(i, 11, 1),pos_eci(i, 11, 2),pos_eci(i, 11, 3),gs_eci(i, 2, 1)-pos_eci(i, 11, 1),gs_eci(i, 2, 2)-pos_eci(i, 11, 2),gs_eci(i, 2, 3)-pos_eci(i, 11, 3),'r');
        end
        visT2_12 =  visible_time{2,12};
        if any(abs(visT2_12-i)<0.0001)
            visL2_12 = quiver3(pos_eci(i, 12, 1),pos_eci(i, 12, 2),pos_eci(i, 12, 3),gs_eci(i, 2, 1)-pos_eci(i, 12, 1),gs_eci(i, 2, 2)-pos_eci(i, 12, 2),gs_eci(i, 2, 3)-pos_eci(i, 12, 3),'r');
        end
        
        % GS 3
        visT3_1 =  visible_time{3,1};
        if any(abs(visT3_1-i)<0.0001)
            visL3_1 = quiver3(pos_eci(i, 1, 1),pos_eci(i, 1, 2),pos_eci(i, 1, 3),gs_eci(i, 3, 1)-pos_eci(i, 1, 1),gs_eci(i, 3, 2)-pos_eci(i, 1, 2),gs_eci(i, 3, 3)-pos_eci(i, 1, 3),'r');
        end
        visT3_2 =  visible_time{3,2};
        if any(abs(visT3_2-i)<0.0001)
            visL3_2 = quiver3(pos_eci(i, 2, 1),pos_eci(i, 2, 2),pos_eci(i, 2, 3),gs_eci(i, 3, 1)-pos_eci(i, 2, 1),gs_eci(i, 3, 2)-pos_eci(i, 2, 2),gs_eci(i, 3, 3)-pos_eci(i, 2, 3),'r');
        end
        visT3_3 =  visible_time{3,3};
        if any(abs(visT3_3-i)<0.0001)
            visL3_3 = quiver3(pos_eci(i, 3, 1),pos_eci(i, 3, 2),pos_eci(i, 3, 3),gs_eci(i, 3, 1)-pos_eci(i, 3, 1),gs_eci(i, 3, 2)-pos_eci(i, 3, 2),gs_eci(i, 3, 3)-pos_eci(i, 3, 3),'r');
        end
        visT3_4 =  visible_time{3,4};
        if any(abs(visT3_4-i)<0.0001)
            visL3_4 = quiver3(pos_eci(i, 4, 1),pos_eci(i, 4, 2),pos_eci(i, 4, 3),gs_eci(i, 3, 1)-pos_eci(i, 4, 1),gs_eci(i, 3, 2)-pos_eci(i, 4, 2),gs_eci(i, 3, 3)-pos_eci(i, 4, 3),'r');
        end
        visT3_5 =  visible_time{3,5};
        if any(abs(visT3_5-i)<0.0001)
            visL3_5 = quiver3(pos_eci(i, 5, 1),pos_eci(i, 5, 2),pos_eci(i, 5, 3),gs_eci(i, 3, 1)-pos_eci(i, 5, 1),gs_eci(i, 3, 2)-pos_eci(i, 5, 2),gs_eci(i, 3, 3)-pos_eci(i, 5, 3),'r');
        end
        visT3_6 =  visible_time{3,6};
        if any(abs(visT3_6-i)<0.0001)
            visL3_6 = quiver3(pos_eci(i, 6, 1),pos_eci(i, 6, 2),pos_eci(i, 6, 3),gs_eci(i, 3, 1)-pos_eci(i, 6, 1),gs_eci(i, 3, 2)-pos_eci(i, 6, 2),gs_eci(i, 3, 3)-pos_eci(i, 6, 3),'r');
        end
        visT3_7 =  visible_time{3,7};
        if any(abs(visT3_7-i)<0.0001)
            visL3_7 = quiver3(pos_eci(i, 7, 1),pos_eci(i, 7, 2),pos_eci(i, 7, 3),gs_eci(i, 3, 1)-pos_eci(i, 7, 1),gs_eci(i, 3, 2)-pos_eci(i, 7, 2),gs_eci(i, 3, 3)-pos_eci(i, 7, 3),'r');
        end
        visT3_8 =  visible_time{3,8};
        if any(abs(visT3_8-i)<0.0001)
            visL3_8 = quiver3(pos_eci(i, 8, 1),pos_eci(i, 8, 2),pos_eci(i, 8, 3),gs_eci(i, 3, 1)-pos_eci(i, 8, 1),gs_eci(i, 3, 2)-pos_eci(i, 8, 2),gs_eci(i, 3, 3)-pos_eci(i, 8, 3),'r');
        end
        visT3_9 =  visible_time{3,9};
        if any(abs(visT3_9-i)<0.0001)
            visL3_9 = quiver3(pos_eci(i, 9, 1),pos_eci(i, 9, 2),pos_eci(i, 9, 3),gs_eci(i, 3, 1)-pos_eci(i, 9, 1),gs_eci(i, 3, 2)-pos_eci(i, 9, 2),gs_eci(i, 3, 3)-pos_eci(i, 9, 3),'r');
        end
        visT3_10 =  visible_time{3,10};
        if any(abs(visT3_10-i)<0.0001)
            visL3_10 = quiver3(pos_eci(i, 10, 1),pos_eci(i, 10, 2),pos_eci(i, 10, 3),gs_eci(i, 3, 1)-pos_eci(i, 10, 1),gs_eci(i, 3, 2)-pos_eci(i, 10, 2),gs_eci(i, 3, 3)-pos_eci(i, 10, 3),'r');
        end
        visT3_11 =  visible_time{3,11};
        if any(abs(visT3_11-i)<0.0001)
            visL3_11 = quiver3(pos_eci(i, 11, 1),pos_eci(i, 11, 2),pos_eci(i, 11, 3),gs_eci(i, 3, 1)-pos_eci(i, 11, 1),gs_eci(i, 3, 2)-pos_eci(i, 11, 2),gs_eci(i, 3, 3)-pos_eci(i, 11, 3),'r');
        end
        visT3_12 =  visible_time{3,12};
        if any(abs(visT3_12-i)<0.0001)
            visL3_12 = quiver3(pos_eci(i, 12, 1),pos_eci(i, 12, 2),pos_eci(i, 12, 3),gs_eci(i, 3, 1)-pos_eci(i, 12, 1),gs_eci(i, 3, 2)-pos_eci(i, 12, 2),gs_eci(i, 3, 3)-pos_eci(i, 12, 3),'r');
        end
        
        % GS 4
        visT4_1 =  visible_time{4,1};
        if any(abs(visT4_1-i)<0.0001)
            visL4_1 = quiver3(pos_eci(i, 1, 1),pos_eci(i, 1, 2),pos_eci(i, 1, 3),gs_eci(i, 4, 1)-pos_eci(i, 1, 1),gs_eci(i, 4, 2)-pos_eci(i, 1, 2),gs_eci(i, 4, 3)-pos_eci(i, 1, 3),'r');
        end
        visT4_2 =  visible_time{4,2};
        if any(abs(visT4_2-i)<0.0001)
            visL4_2 = quiver3(pos_eci(i, 2, 1),pos_eci(i, 2, 2),pos_eci(i, 2, 3),gs_eci(i, 4, 1)-pos_eci(i, 2, 1),gs_eci(i, 4, 2)-pos_eci(i, 2, 2),gs_eci(i, 4, 3)-pos_eci(i, 2, 3),'r');
        end
        visT4_3 =  visible_time{4,3};
        if any(abs(visT4_3-i)<0.0001)
            visL4_3 = quiver3(pos_eci(i, 3, 1),pos_eci(i, 3, 2),pos_eci(i, 3, 3),gs_eci(i, 4, 1)-pos_eci(i, 3, 1),gs_eci(i, 4, 2)-pos_eci(i, 3, 2),gs_eci(i, 4, 3)-pos_eci(i, 3, 3),'r');
        end
        visT4_4 =  visible_time{4,4};
        if any(abs(visT4_4-i)<0.0001)
            visL4_4 = quiver3(pos_eci(i, 4, 1),pos_eci(i, 4, 2),pos_eci(i, 4, 3),gs_eci(i, 4, 1)-pos_eci(i, 4, 1),gs_eci(i, 4, 2)-pos_eci(i, 4, 2),gs_eci(i, 4, 3)-pos_eci(i, 4, 3),'r');
        end
        visT4_5 =  visible_time{4,5};
        if any(abs(visT4_5-i)<0.0001)
            visL4_5 = quiver3(pos_eci(i, 5, 1),pos_eci(i, 5, 2),pos_eci(i, 5, 3),gs_eci(i, 4, 1)-pos_eci(i, 5, 1),gs_eci(i, 4, 2)-pos_eci(i, 5, 2),gs_eci(i, 4, 3)-pos_eci(i, 5, 3),'r');
        end
        visT4_6 =  visible_time{4,6};
        if any(abs(visT4_6-i)<0.0001)
            visL4_6 = quiver3(pos_eci(i, 6, 1),pos_eci(i, 6, 2),pos_eci(i, 6, 3),gs_eci(i, 4, 1)-pos_eci(i, 6, 1),gs_eci(i, 4, 2)-pos_eci(i, 6, 2),gs_eci(i, 4, 3)-pos_eci(i, 6, 3),'r');
        end
        visT4_7 =  visible_time{4,7};
        if any(abs(visT4_7-i)<0.0001)
            visL4_7 = quiver3(pos_eci(i, 7, 1),pos_eci(i, 7, 2),pos_eci(i, 7, 3),gs_eci(i, 4, 1)-pos_eci(i, 7, 1),gs_eci(i, 4, 2)-pos_eci(i, 7, 2),gs_eci(i, 4, 3)-pos_eci(i, 7, 3),'r');
        end
        visT4_8 =  visible_time{4,8};
        if any(abs(visT4_8-i)<0.0001)
            visL4_8 = quiver3(pos_eci(i, 8, 1),pos_eci(i, 8, 2),pos_eci(i, 8, 3),gs_eci(i, 4, 1)-pos_eci(i, 8, 1),gs_eci(i, 4, 2)-pos_eci(i, 8, 2),gs_eci(i, 4, 3)-pos_eci(i, 8, 3),'r');
        end
        visT4_9 =  visible_time{4,9};
        if any(abs(visT4_9-i)<0.0001)
            visL4_9 = quiver3(pos_eci(i, 9, 1),pos_eci(i, 9, 2),pos_eci(i, 9, 3),gs_eci(i, 4, 1)-pos_eci(i, 9, 1),gs_eci(i, 4, 2)-pos_eci(i, 9, 2),gs_eci(i, 4, 3)-pos_eci(i, 9, 3),'r');
        end
        visT4_10 =  visible_time{4,10};
        if any(abs(visT4_10-i)<0.0001)
            visL4_10 = quiver3(pos_eci(i, 10, 1),pos_eci(i, 10, 2),pos_eci(i, 10, 3),gs_eci(i, 4, 1)-pos_eci(i, 10, 1),gs_eci(i, 4, 2)-pos_eci(i, 10, 2),gs_eci(i, 4, 3)-pos_eci(i, 10, 3),'r');
        end
        visT4_11 =  visible_time{4,11};
        if any(abs(visT4_11-i)<0.0001)
            visL4_11 = quiver3(pos_eci(i, 11, 1),pos_eci(i, 11, 2),pos_eci(i, 11, 3),gs_eci(i, 4, 1)-pos_eci(i, 11, 1),gs_eci(i, 4, 2)-pos_eci(i, 11, 2),gs_eci(i, 4, 3)-pos_eci(i, 11, 3),'r');
        end
        visT4_12 =  visible_time{4,12};
        if any(abs(visT4_12-i)<0.0001)
            visL4_12 = quiver3(pos_eci(i, 12, 1),pos_eci(i, 12, 2),pos_eci(i, 12, 3),gs_eci(i, 4, 1)-pos_eci(i, 12, 1),gs_eci(i, 4, 2)-pos_eci(i, 12, 2),gs_eci(i, 4, 3)-pos_eci(i, 12, 3),'r');
        end
        
        % GS 5
        visT5_1 =  visible_time{5,1};
        if any(abs(visT5_1-i)<0.0001)
            visL5_1 = quiver3(pos_eci(i, 1, 1),pos_eci(i, 1, 2),pos_eci(i, 1, 3),gs_eci(i, 5, 1)-pos_eci(i, 1, 1),gs_eci(i, 5, 2)-pos_eci(i, 1, 2),gs_eci(i, 5, 3)-pos_eci(i, 1, 3),'r');
        end
        visT5_2 =  visible_time{5,2};
        if any(abs(visT5_2-i)<0.0001)
            visL5_2 = quiver3(pos_eci(i, 2, 1),pos_eci(i, 2, 2),pos_eci(i, 2, 3),gs_eci(i, 5, 1)-pos_eci(i, 2, 1),gs_eci(i, 5, 2)-pos_eci(i, 2, 2),gs_eci(i, 5, 3)-pos_eci(i, 2, 3),'r');
        end
        visT5_3 =  visible_time{5,3};
        if any(abs(visT5_3-i)<0.0001)
            visL5_3 = quiver3(pos_eci(i, 3, 1),pos_eci(i, 3, 2),pos_eci(i, 3, 3),gs_eci(i, 5, 1)-pos_eci(i, 3, 1),gs_eci(i, 5, 2)-pos_eci(i, 3, 2),gs_eci(i, 5, 3)-pos_eci(i, 3, 3),'r');
        end
        visT5_4 =  visible_time{5,4};
        if any(abs(visT5_4-i)<0.0001)
            visL5_4 = quiver3(pos_eci(i, 4, 1),pos_eci(i, 4, 2),pos_eci(i, 4, 3),gs_eci(i, 5, 1)-pos_eci(i, 4, 1),gs_eci(i, 5, 2)-pos_eci(i, 4, 2),gs_eci(i, 5, 3)-pos_eci(i, 4, 3),'r');
        end
        visT5_5 =  visible_time{5,5};
        if any(abs(visT5_5-i)<0.0001)
            visL5_5 = quiver3(pos_eci(i, 5, 1),pos_eci(i, 5, 2),pos_eci(i, 5, 3),gs_eci(i, 5, 1)-pos_eci(i, 5, 1),gs_eci(i, 5, 2)-pos_eci(i, 5, 2),gs_eci(i, 5, 3)-pos_eci(i, 5, 3),'r');
        end
        visT5_6 =  visible_time{5,6};
        if any(abs(visT5_6-i)<0.0001)
            visL5_6 = quiver3(pos_eci(i, 6, 1),pos_eci(i, 6, 2),pos_eci(i, 6, 3),gs_eci(i, 5, 1)-pos_eci(i, 6, 1),gs_eci(i, 5, 2)-pos_eci(i, 6, 2),gs_eci(i, 5, 3)-pos_eci(i, 6, 3),'r');
        end
        visT5_7 =  visible_time{5,7};
        if any(abs(visT5_7-i)<0.0001)
            visL5_7 = quiver3(pos_eci(i, 7, 1),pos_eci(i, 7, 2),pos_eci(i, 7, 3),gs_eci(i, 5, 1)-pos_eci(i, 7, 1),gs_eci(i, 5, 2)-pos_eci(i, 7, 2),gs_eci(i, 5, 3)-pos_eci(i, 7, 3),'r');
        end
        visT5_8 =  visible_time{5,8};
        if any(abs(visT5_8-i)<0.0001)
            visL5_8 = quiver3(pos_eci(i, 8, 1),pos_eci(i, 8, 2),pos_eci(i, 8, 3),gs_eci(i, 5, 1)-pos_eci(i, 8, 1),gs_eci(i, 5, 2)-pos_eci(i, 8, 2),gs_eci(i, 5, 3)-pos_eci(i, 8, 3),'r');
        end
        visT5_9 =  visible_time{5,9};
        if any(abs(visT5_9-i)<0.0001)
            visL5_9 = quiver3(pos_eci(i, 9, 1),pos_eci(i, 9, 2),pos_eci(i, 9, 3),gs_eci(i, 5, 1)-pos_eci(i, 9, 1),gs_eci(i, 5, 2)-pos_eci(i, 9, 2),gs_eci(i, 5, 3)-pos_eci(i, 9, 3),'r');
        end
        visT5_10 =  visible_time{5,10};
        if any(abs(visT5_10-i)<0.0001)
            visL5_10 = quiver3(pos_eci(i, 10, 1),pos_eci(i, 10, 2),pos_eci(i, 10, 3),gs_eci(i, 5, 1)-pos_eci(i, 10, 1),gs_eci(i, 5, 2)-pos_eci(i, 10, 2),gs_eci(i, 5, 3)-pos_eci(i, 10, 3),'r');
        end
        visT5_11 =  visible_time{5,11};
        if any(abs(visT5_11-i)<0.0001)
            visL5_11 = quiver3(pos_eci(i, 11, 1),pos_eci(i, 11, 2),pos_eci(i, 11, 3),gs_eci(i, 5, 1)-pos_eci(i, 11, 1),gs_eci(i, 5, 2)-pos_eci(i, 11, 2),gs_eci(i, 5, 3)-pos_eci(i, 11, 3),'r');
        end
        visT5_12 =  visible_time{5,12};
        if any(abs(visT5_12-i)<0.0001)
            visL5_12 = quiver3(pos_eci(i, 12, 1),pos_eci(i, 12, 2),pos_eci(i, 12, 3),gs_eci(i, 5, 1)-pos_eci(i, 12, 1),gs_eci(i, 5, 2)-pos_eci(i, 12, 2),gs_eci(i, 5, 3)-pos_eci(i, 12, 3),'r');
        end
        
    else
        refreshdata(plt.satsECEF, 'caller');
        
        % GS 1
        visT1_1 =  visible_time{1,1};
        if any(abs(visT1_1-i)<0.0001)
            visL1_1 = quiver3(pos_ecef(i, 1, 1),pos_ecef(i, 1, 2),pos_ecef(i, 1, 3),gs_ecef(1, 1)-pos_ecef(i, 1, 1),gs_ecef(2, 1)-pos_ecef(i, 1, 2),gs_ecef( 3, 1)-pos_ecef(i, 1, 3),'r');
        end
        visT1_2 =  visible_time{1,2};
        if any(abs(visT1_2-i)<0.0001)
            visL1_2 = quiver3(pos_ecef(i, 2, 1),pos_ecef(i, 2, 2),pos_ecef(i, 2, 3),gs_ecef(1, 1)-pos_ecef(i, 2, 1),gs_ecef(2, 1)-pos_ecef(i, 2, 2),gs_ecef( 3, 1)-pos_ecef(i, 2, 3),'r');
        end
        visT1_3 =  visible_time{1,3};
        if any(abs(visT1_3-i)<0.0001)
            visL1_3 = quiver3(pos_ecef(i, 3, 1),pos_ecef(i, 3, 2),pos_ecef(i, 3, 3),gs_ecef(1, 1)-pos_ecef(i, 3, 1),gs_ecef(2, 1)-pos_ecef(i, 3, 2),gs_ecef( 3, 1)-pos_ecef(i, 3, 3),'r');
        end
        visT1_4 =  visible_time{1,4};
        if any(abs(visT1_4-i)<0.0001)
            visL1_4 = quiver3(pos_ecef(i, 4, 1),pos_ecef(i, 4, 2),pos_ecef(i, 4, 3),gs_ecef(1, 1)-pos_ecef(i, 4, 1),gs_ecef(2, 1)-pos_ecef(i, 4, 2),gs_ecef( 3, 1)-pos_ecef(i, 4, 3),'r');
        end
        visT1_5 =  visible_time{1,5};
        if any(abs(visT1_5-i)<0.0001)
            visL1_5 = quiver3(pos_ecef(i, 5, 1),pos_ecef(i, 5, 2),pos_ecef(i, 5, 3),gs_ecef(1, 1)-pos_ecef(i, 5, 1),gs_ecef(2, 1)-pos_ecef(i, 5, 2),gs_ecef( 3, 1)-pos_ecef(i, 5, 3),'r');
        end
        visT1_6 =  visible_time{1,6};
        if any(abs(visT1_6-i)<0.0001)
            visL1_6 = quiver3(pos_ecef(i, 6, 1),pos_ecef(i, 6, 2),pos_ecef(i, 6, 3),gs_ecef(1, 1)-pos_ecef(i, 6, 1),gs_ecef(2, 1)-pos_ecef(i, 6, 2),gs_ecef( 3, 1)-pos_ecef(i, 6, 3),'r');
        end
        visT1_7 =  visible_time{1,7};
        if any(abs(visT1_7-i)<0.0001)
            visL1_7 = quiver3(pos_ecef(i, 7, 1),pos_ecef(i, 7, 2),pos_ecef(i, 7, 3),gs_ecef(1, 1)-pos_ecef(i, 7, 1),gs_ecef(2, 1)-pos_ecef(i, 7, 2),gs_ecef( 3, 1)-pos_ecef(i, 7, 3),'r');
        end
        visT1_8 =  visible_time{1,8};
        if any(abs(visT1_8-i)<0.0001)
            visL1_8 = quiver3(pos_ecef(i, 8, 1),pos_ecef(i, 8, 2),pos_ecef(i, 8, 3),gs_ecef(1, 1)-pos_ecef(i, 8, 1),gs_ecef(2, 1)-pos_ecef(i, 8, 2),gs_ecef( 3, 1)-pos_ecef(i, 8, 3),'r');
        end
        visT1_9 =  visible_time{1,9};
        if any(abs(visT1_9-i)<0.0001)
            visL1_9 = quiver3(pos_ecef(i, 9, 1),pos_ecef(i, 9, 2),pos_ecef(i, 9, 3),gs_ecef(1, 1)-pos_ecef(i, 9, 1),gs_ecef(2, 1)-pos_ecef(i, 9, 2),gs_ecef( 3, 1)-pos_ecef(i, 9, 3),'r');
        end
        visT1_10 =  visible_time{1,10};
        if any(abs(visT1_10-i)<0.0001)
            visL1_10 = quiver3(pos_ecef(i, 10, 1),pos_ecef(i, 10, 2),pos_ecef(i, 10, 3),gs_ecef(1, 1)-pos_ecef(i, 10, 1),gs_ecef(2, 1)-pos_ecef(i, 10, 2),gs_ecef( 3, 1)-pos_ecef(i, 10, 3),'r');
        end
        visT1_11 =  visible_time{1,11};
        if any(abs(visT1_11-i)<0.0001)
            visL1_11 = quiver3(pos_ecef(i, 11, 1),pos_ecef(i, 11, 2),pos_ecef(i, 11, 3),gs_ecef(1, 1)-pos_ecef(i, 11, 1),gs_ecef(2, 1)-pos_ecef(i, 11, 2),gs_ecef( 3, 1)-pos_ecef(i, 11, 3),'r');
        end
        visT1_12 =  visible_time{1,12};
        if any(abs(visT1_12-i)<0.0001)
            visL1_12 = quiver3(pos_ecef(i, 12, 1),pos_ecef(i, 12, 2),pos_ecef(i, 12, 3),gs_ecef(1, 1)-pos_ecef(i, 12, 1),gs_ecef(2, 1)-pos_ecef(i, 12, 2),gs_ecef( 3, 1)-pos_ecef(i, 12, 3),'r');
        end
        
        % GS 2
        visT2_1 =  visible_time{2,1};
        if any(abs(visT2_1-i)<0.0001)
            visL2_1 = quiver3(pos_ecef(i, 1, 1),pos_ecef(i, 1, 2),pos_ecef(i, 1, 3),gs_ecef(1, 2)-pos_ecef(i, 1, 1),gs_ecef(2, 2)-pos_ecef(i, 1, 2),gs_ecef( 3, 2)-pos_ecef(i, 1, 3),'r');
        end
        visT2_2 =  visible_time{2,2};
        if any(abs(visT2_2-i)<0.0001)
            visL2_2 = quiver3(pos_ecef(i, 2, 1),pos_ecef(i, 2, 2),pos_ecef(i, 2, 3),gs_ecef(1, 2)-pos_ecef(i, 2, 1),gs_ecef(2, 2)-pos_ecef(i, 2, 2),gs_ecef( 3, 2)-pos_ecef(i, 2, 3),'r');
        end
        visT2_3 =  visible_time{2,3};
        if any(abs(visT2_3-i)<0.0001)
            visL2_3 = quiver3(pos_ecef(i, 3, 1),pos_ecef(i, 3, 2),pos_ecef(i, 3, 3),gs_ecef(1, 2)-pos_ecef(i, 3, 1),gs_ecef(2, 2)-pos_ecef(i, 3, 2),gs_ecef( 3, 2)-pos_ecef(i, 3, 3),'r');
        end
        visT2_4 =  visible_time{2,4};
        if any(abs(visT2_4-i)<0.0001)
            visL2_4 = quiver3(pos_ecef(i, 4, 1),pos_ecef(i, 4, 2),pos_ecef(i, 4, 3),gs_ecef(1, 2)-pos_ecef(i, 4, 1),gs_ecef(2, 2)-pos_ecef(i, 4, 2),gs_ecef( 3, 2)-pos_ecef(i, 4, 3),'r');
        end
        visT2_5 =  visible_time{2,5};
        if any(abs(visT2_5-i)<0.0001)
            visL2_5 = quiver3(pos_ecef(i, 5, 1),pos_ecef(i, 5, 2),pos_ecef(i, 5, 3),gs_ecef(1, 2)-pos_ecef(i, 5, 1),gs_ecef(2, 2)-pos_ecef(i, 5, 2),gs_ecef( 3, 2)-pos_ecef(i, 5, 3),'r');
        end
        visT2_6 =  visible_time{2,6};
        if any(abs(visT2_6-i)<0.0001)
            visL2_6 = quiver3(pos_ecef(i, 6, 1),pos_ecef(i, 6, 2),pos_ecef(i, 6, 3),gs_ecef(1, 2)-pos_ecef(i, 6, 1),gs_ecef(2, 2)-pos_ecef(i, 6, 2),gs_ecef( 3, 2)-pos_ecef(i, 6, 3),'r');
        end
        visT2_7 =  visible_time{2,7};
        if any(abs(visT2_7-i)<0.0001)
            visL2_7 = quiver3(pos_ecef(i, 7, 1),pos_ecef(i, 7, 2),pos_ecef(i, 7, 3),gs_ecef(1, 2)-pos_ecef(i, 7, 1),gs_ecef(2, 2)-pos_ecef(i, 7, 2),gs_ecef( 3, 2)-pos_ecef(i, 7, 3),'r');
        end
        visT2_8 =  visible_time{2,8};
        if any(abs(visT2_8-i)<0.0001)
            visL2_8 = quiver3(pos_ecef(i, 8, 1),pos_ecef(i, 8, 2),pos_ecef(i, 8, 3),gs_ecef(1, 2)-pos_ecef(i, 8, 1),gs_ecef(2, 2)-pos_ecef(i, 8, 2),gs_ecef( 3, 2)-pos_ecef(i, 8, 3),'r');
        end
        visT2_9 =  visible_time{2,9};
        if any(abs(visT2_9-i)<0.0001)
            visL2_9 = quiver3(pos_ecef(i, 9, 1),pos_ecef(i, 9, 2),pos_ecef(i, 9, 3),gs_ecef(1, 2)-pos_ecef(i, 9, 1),gs_ecef(2, 2)-pos_ecef(i, 9, 2),gs_ecef( 3, 2)-pos_ecef(i, 9, 3),'r');
        end
        visT2_10 =  visible_time{2,10};
        if any(abs(visT2_10-i)<0.0001)
            visL2_10 = quiver3(pos_ecef(i, 10, 1),pos_ecef(i, 10, 2),pos_ecef(i, 10, 3),gs_ecef(1, 2)-pos_ecef(i, 10, 1),gs_ecef(2, 2)-pos_ecef(i, 10, 2),gs_ecef( 3, 2)-pos_ecef(i, 10, 3),'r');
        end
        visT2_11 =  visible_time{2,11};
        if any(abs(visT2_11-i)<0.0001)
            visL2_11 = quiver3(pos_ecef(i, 11, 1),pos_ecef(i, 11, 2),pos_ecef(i, 11, 3),gs_ecef(1, 2)-pos_ecef(i, 11, 1),gs_ecef(2, 2)-pos_ecef(i, 11, 2),gs_ecef( 3, 2)-pos_ecef(i, 11, 3),'r');
        end
        visT2_12 =  visible_time{2,12};
        if any(abs(visT2_12-i)<0.0001)
            visL2_12 = quiver3(pos_ecef(i, 12, 1),pos_ecef(i, 12, 2),pos_ecef(i, 12, 3),gs_ecef(1, 2)-pos_ecef(i, 12, 1),gs_ecef(2, 2)-pos_ecef(i, 12, 2),gs_ecef( 3, 2)-pos_ecef(i, 12, 3),'r');
        end
        
        % GS 3
        visT3_1 =  visible_time{3,1};
        if any(abs(visT3_1-i)<0.0001)
            visL3_1 = quiver3(pos_ecef(i, 1, 1),pos_ecef(i, 1, 2),pos_ecef(i, 1, 3),gs_ecef(1, 3)-pos_ecef(i, 1, 1),gs_ecef(2, 3)-pos_ecef(i, 1, 2),gs_ecef( 3, 3)-pos_ecef(i, 1, 3),'r');
        end
        visT3_2 =  visible_time{3,2};
        if any(abs(visT3_2-i)<0.0001)
            visL3_2 = quiver3(pos_ecef(i, 2, 1),pos_ecef(i, 2, 2),pos_ecef(i, 2, 3),gs_ecef(1, 3)-pos_ecef(i, 2, 1),gs_ecef(2, 3)-pos_ecef(i, 2, 2),gs_ecef( 3, 3)-pos_ecef(i, 2, 3),'r');
        end
        visT3_3 =  visible_time{3,3};
        if any(abs(visT3_3-i)<0.0001)
            visL3_3 = quiver3(pos_ecef(i, 3, 1),pos_ecef(i, 3, 2),pos_ecef(i, 3, 3),gs_ecef(1, 3)-pos_ecef(i, 3, 1),gs_ecef(2, 3)-pos_ecef(i, 3, 2),gs_ecef( 3, 3)-pos_ecef(i, 3, 3),'r');
        end
        visT3_4 =  visible_time{3,4};
        if any(abs(visT3_4-i)<0.0001)
            visL3_4 = quiver3(pos_ecef(i, 4, 1),pos_ecef(i, 4, 2),pos_ecef(i, 4, 3),gs_ecef(1, 3)-pos_ecef(i, 4, 1),gs_ecef(2, 3)-pos_ecef(i, 4, 2),gs_ecef( 3, 3)-pos_ecef(i, 4, 3),'r');
        end
        visT3_5 =  visible_time{3,5};
        if any(abs(visT3_5-i)<0.0001)
            visL3_5 = quiver3(pos_ecef(i, 5, 1),pos_ecef(i, 5, 2),pos_ecef(i, 5, 3),gs_ecef(1, 3)-pos_ecef(i, 5, 1),gs_ecef(2, 3)-pos_ecef(i, 5, 2),gs_ecef( 3, 3)-pos_ecef(i, 5, 3),'r');
        end
        visT3_6 =  visible_time{3,6};
        if any(abs(visT3_6-i)<0.0001)
            visL3_6 = quiver3(pos_ecef(i, 6, 1),pos_ecef(i, 6, 2),pos_ecef(i, 6, 3),gs_ecef(1, 3)-pos_ecef(i, 6, 1),gs_ecef(2, 3)-pos_ecef(i, 6, 2),gs_ecef( 3, 3)-pos_ecef(i, 6, 3),'r');
        end
        visT3_7 =  visible_time{3,7};
        if any(abs(visT3_7-i)<0.0001)
            visL3_7 = quiver3(pos_ecef(i, 7, 1),pos_ecef(i, 7, 2),pos_ecef(i, 7, 3),gs_ecef(1, 3)-pos_ecef(i, 7, 1),gs_ecef(2, 3)-pos_ecef(i, 7, 2),gs_ecef( 3, 3)-pos_ecef(i, 7, 3),'r');
        end
        visT3_8 =  visible_time{3,8};
        if any(abs(visT3_8-i)<0.0001)
            visL3_8 = quiver3(pos_ecef(i, 8, 1),pos_ecef(i, 8, 2),pos_ecef(i, 8, 3),gs_ecef(1, 3)-pos_ecef(i, 8, 1),gs_ecef(2, 3)-pos_ecef(i, 8, 2),gs_ecef( 3, 3)-pos_ecef(i, 8, 3),'r');
        end
        visT3_9 =  visible_time{3,9};
        if any(abs(visT3_9-i)<0.0001)
            visL3_9 = quiver3(pos_ecef(i, 9, 1),pos_ecef(i, 9, 2),pos_ecef(i, 9, 3),gs_ecef(1, 3)-pos_ecef(i, 9, 1),gs_ecef(2, 3)-pos_ecef(i, 9, 2),gs_ecef( 3, 3)-pos_ecef(i, 9, 3),'r');
        end
        visT3_10 =  visible_time{3,10};
        if any(abs(visT3_10-i)<0.0001)
            visL3_10 = quiver3(pos_ecef(i, 10, 1),pos_ecef(i, 10, 2),pos_ecef(i, 10, 3),gs_ecef(1, 3)-pos_ecef(i, 10, 1),gs_ecef(2, 3)-pos_ecef(i, 10, 2),gs_ecef( 3, 3)-pos_ecef(i, 10, 3),'r');
        end
        visT3_11 =  visible_time{3,11};
        if any(abs(visT3_11-i)<0.0001)
            visL3_11 = quiver3(pos_ecef(i, 11, 1),pos_ecef(i, 11, 2),pos_ecef(i, 11, 3),gs_ecef(1, 3)-pos_ecef(i, 11, 1),gs_ecef(2, 3)-pos_ecef(i, 11, 2),gs_ecef( 3, 3)-pos_ecef(i, 11, 3),'r');
        end
        visT3_12 =  visible_time{3,12};
        if any(abs(visT3_12-i)<0.0001)
            visL3_12 = quiver3(pos_ecef(i, 12, 1),pos_ecef(i, 12, 2),pos_ecef(i, 12, 3),gs_ecef(1, 3)-pos_ecef(i, 12, 1),gs_ecef(2, 3)-pos_ecef(i, 12, 2),gs_ecef( 3, 3)-pos_ecef(i, 12, 3),'r');
        end
        
        % GS 4
        visT4_1 =  visible_time{4,1};
        if any(abs(visT4_1-i)<0.0001)
            visL4_1 = quiver3(pos_ecef(i, 1, 1),pos_ecef(i, 1, 2),pos_ecef(i, 1, 3),gs_ecef(1, 4)-pos_ecef(i, 1, 1),gs_ecef(2, 4)-pos_ecef(i, 1, 2),gs_ecef( 3, 4)-pos_ecef(i, 1, 3),'r');
        end
        visT4_2 =  visible_time{4,2};
        if any(abs(visT4_2-i)<0.0001)
            visL4_2 = quiver3(pos_ecef(i, 2, 1),pos_ecef(i, 2, 2),pos_ecef(i, 2, 3),gs_ecef(1, 4)-pos_ecef(i, 2, 1),gs_ecef(2, 4)-pos_ecef(i, 2, 2),gs_ecef( 3, 4)-pos_ecef(i, 2, 3),'r');
        end
        visT4_3 =  visible_time{4,3};
        if any(abs(visT4_3-i)<0.0001)
            visL4_3 = quiver3(pos_ecef(i, 3, 1),pos_ecef(i, 3, 2),pos_ecef(i, 3, 3),gs_ecef(1, 4)-pos_ecef(i, 3, 1),gs_ecef(2, 4)-pos_ecef(i, 3, 2),gs_ecef( 3, 4)-pos_ecef(i, 3, 3),'r');
        end
        visT4_4 =  visible_time{4,4};
        if any(abs(visT4_4-i)<0.0001)
            visL4_4 = quiver3(pos_ecef(i, 4, 1),pos_ecef(i, 4, 2),pos_ecef(i, 4, 3),gs_ecef(1, 4)-pos_ecef(i, 4, 1),gs_ecef(2, 4)-pos_ecef(i, 4, 2),gs_ecef( 3, 4)-pos_ecef(i, 4, 3),'r');
        end
        visT4_5 =  visible_time{4,5};
        if any(abs(visT4_5-i)<0.0001)
            visL4_5 = quiver3(pos_ecef(i, 5, 1),pos_ecef(i, 5, 2),pos_ecef(i, 5, 3),gs_ecef(1, 4)-pos_ecef(i, 5, 1),gs_ecef(2, 4)-pos_ecef(i, 5, 2),gs_ecef( 3, 4)-pos_ecef(i, 5, 3),'r');
        end
        visT4_6 =  visible_time{4,6};
        if any(abs(visT4_6-i)<0.0001)
            visL4_6 = quiver3(pos_ecef(i, 6, 1),pos_ecef(i, 6, 2),pos_ecef(i, 6, 3),gs_ecef(1, 4)-pos_ecef(i, 6, 1),gs_ecef(2, 4)-pos_ecef(i, 6, 2),gs_ecef( 3, 4)-pos_ecef(i, 6, 3),'r');
        end
        visT4_7 =  visible_time{4,7};
        if any(abs(visT4_7-i)<0.0001)
            visL4_7 = quiver3(pos_ecef(i, 7, 1),pos_ecef(i, 7, 2),pos_ecef(i, 7, 3),gs_ecef(1, 4)-pos_ecef(i, 7, 1),gs_ecef(2, 4)-pos_ecef(i, 7, 2),gs_ecef( 3, 4)-pos_ecef(i, 7, 3),'r');
        end
        visT4_8 =  visible_time{4,8};
        if any(abs(visT4_8-i)<0.0001)
            visL4_8 = quiver3(pos_ecef(i, 8, 1),pos_ecef(i, 8, 2),pos_ecef(i, 8, 3),gs_ecef(1, 4)-pos_ecef(i, 8, 1),gs_ecef(2, 4)-pos_ecef(i, 8, 2),gs_ecef( 3, 4)-pos_ecef(i, 8, 3),'r');
        end
        visT4_9 =  visible_time{4,9};
        if any(abs(visT4_9-i)<0.0001)
            visL4_9 = quiver3(pos_ecef(i, 9, 1),pos_ecef(i, 9, 2),pos_ecef(i, 9, 3),gs_ecef(1, 4)-pos_ecef(i, 9, 1),gs_ecef(2, 4)-pos_ecef(i, 9, 2),gs_ecef( 3, 4)-pos_ecef(i, 9, 3),'r');
        end
        visT4_10 =  visible_time{4,10};
        if any(abs(visT4_10-i)<0.0001)
            visL4_10 = quiver3(pos_ecef(i, 10, 1),pos_ecef(i, 10, 2),pos_ecef(i, 10, 3),gs_ecef(1, 4)-pos_ecef(i, 10, 1),gs_ecef(2, 4)-pos_ecef(i, 10, 2),gs_ecef( 3, 4)-pos_ecef(i, 10, 3),'r');
        end
        visT4_11 =  visible_time{4,11};
        if any(abs(visT4_11-i)<0.0001)
            visL4_11 = quiver3(pos_ecef(i, 11, 1),pos_ecef(i, 11, 2),pos_ecef(i, 11, 3),gs_ecef(1, 4)-pos_ecef(i, 11, 1),gs_ecef(2, 4)-pos_ecef(i, 11, 2),gs_ecef( 3, 4)-pos_ecef(i, 11, 3),'r');
        end
        visT4_12 =  visible_time{4,12};
        if any(abs(visT4_12-i)<0.0001)
            visL4_12 = quiver3(pos_ecef(i, 12, 1),pos_ecef(i, 12, 2),pos_ecef(i, 12, 3),gs_ecef(1, 4)-pos_ecef(i, 12, 1),gs_ecef(2, 4)-pos_ecef(i, 12, 2),gs_ecef( 3, 4)-pos_ecef(i, 12, 3),'r');
        end
        
        % GS 5
        visT5_1 =  visible_time{5,1};
        if any(abs(visT5_1-i)<0.0001)
            visL5_1 = quiver3(pos_ecef(i, 1, 1),pos_ecef(i, 1, 2),pos_ecef(i, 1, 3),gs_ecef(1, 5)-pos_ecef(i, 1, 1),gs_ecef(2, 5)-pos_ecef(i, 1, 2),gs_ecef( 3, 5)-pos_ecef(i, 1, 3),'r');
        end
        visT5_2 =  visible_time{5,2};
        if any(abs(visT5_2-i)<0.0001)
            visL5_2 = quiver3(pos_ecef(i, 2, 1),pos_ecef(i, 2, 2),pos_ecef(i, 2, 3),gs_ecef(1, 5)-pos_ecef(i, 2, 1),gs_ecef(2, 5)-pos_ecef(i, 2, 2),gs_ecef( 3, 5)-pos_ecef(i, 2, 3),'r');
        end
        visT5_3 =  visible_time{5,3};
        if any(abs(visT5_3-i)<0.0001)
            visL5_3 = quiver3(pos_ecef(i, 3, 1),pos_ecef(i, 3, 2),pos_ecef(i, 3, 3),gs_ecef(1, 5)-pos_ecef(i, 3, 1),gs_ecef(2, 5)-pos_ecef(i, 3, 2),gs_ecef( 3, 5)-pos_ecef(i, 3, 3),'r');
        end
        visT5_4 =  visible_time{5,4};
        if any(abs(visT5_4-i)<0.0001)
            visL5_4 = quiver3(pos_ecef(i, 4, 1),pos_ecef(i, 4, 2),pos_ecef(i, 4, 3),gs_ecef(1, 5)-pos_ecef(i, 4, 1),gs_ecef(2, 5)-pos_ecef(i, 4, 2),gs_ecef( 3, 5)-pos_ecef(i, 4, 3),'r');
        end
        visT5_5 =  visible_time{5,5};
        if any(abs(visT5_5-i)<0.0001)
            visL5_5 = quiver3(pos_ecef(i, 5, 1),pos_ecef(i, 5, 2),pos_ecef(i, 5, 3),gs_ecef(1, 5)-pos_ecef(i, 5, 1),gs_ecef(2, 5)-pos_ecef(i, 5, 2),gs_ecef( 3, 5)-pos_ecef(i, 5, 3),'r');
        end
        visT5_6 =  visible_time{5,6};
        if any(abs(visT5_6-i)<0.0001)
            visL5_6 = quiver3(pos_ecef(i, 6, 1),pos_ecef(i, 6, 2),pos_ecef(i, 6, 3),gs_ecef(1, 5)-pos_ecef(i, 6, 1),gs_ecef(2, 5)-pos_ecef(i, 6, 2),gs_ecef( 3, 5)-pos_ecef(i, 6, 3),'r');
        end
        visT5_7 =  visible_time{5,7};
        if any(abs(visT5_7-i)<0.0001)
            visL5_7 = quiver3(pos_ecef(i, 7, 1),pos_ecef(i, 7, 2),pos_ecef(i, 7, 3),gs_ecef(1, 5)-pos_ecef(i, 7, 1),gs_ecef(2, 5)-pos_ecef(i, 7, 2),gs_ecef( 3, 5)-pos_ecef(i, 7, 3),'r');
        end
        visT5_8 =  visible_time{5,8};
        if any(abs(visT5_8-i)<0.0001)
            visL5_8 = quiver3(pos_ecef(i, 8, 1),pos_ecef(i, 8, 2),pos_ecef(i, 8, 3),gs_ecef(1, 5)-pos_ecef(i, 8, 1),gs_ecef(2, 5)-pos_ecef(i, 8, 2),gs_ecef( 3, 5)-pos_ecef(i, 8, 3),'r');
        end
        visT5_9 =  visible_time{5,9};
        if any(abs(visT5_9-i)<0.0001)
            visL5_9 = quiver3(pos_ecef(i, 9, 1),pos_ecef(i, 9, 2),pos_ecef(i, 9, 3),gs_ecef(1, 5)-pos_ecef(i, 9, 1),gs_ecef(2, 5)-pos_ecef(i, 9, 2),gs_ecef( 3, 5)-pos_ecef(i, 9, 3),'r');
        end
        visT5_10 =  visible_time{5,10};
        if any(abs(visT5_10-i)<0.0001)
            visL5_10 = quiver3(pos_ecef(i, 10, 1),pos_ecef(i, 10, 2),pos_ecef(i, 10, 3),gs_ecef(1, 5)-pos_ecef(i, 10, 1),gs_ecef(2, 5)-pos_ecef(i, 10, 2),gs_ecef( 3, 5)-pos_ecef(i, 10, 3),'r');
        end
        visT5_11 =  visible_time{5,11};
        if any(abs(visT5_11-i)<0.0001)
            visL5_11 = quiver3(pos_ecef(i, 11, 1),pos_ecef(i, 11, 2),pos_ecef(i, 11, 3),gs_ecef(1, 5)-pos_ecef(i, 11, 1),gs_ecef(2, 5)-pos_ecef(i, 11, 2),gs_ecef( 3, 5)-pos_ecef(i, 11, 3),'r');
        end
        visT5_12 =  visible_time{5,12};
        if any(abs(visT5_12-i)<0.0001)
            visL5_12 = quiver3(pos_ecef(i, 12, 1),pos_ecef(i, 12, 2),pos_ecef(i, 12, 3),gs_ecef(1, 5)-pos_ecef(i, 12, 1),gs_ecef(2, 5)-pos_ecef(i, 12, 2),gs_ecef( 3, 5)-pos_ecef(i, 12, 3),'r');
        end
        
    end
    refreshdata(plt.ground_trace, 'caller');
    drawnow();
    
    pause(0.0001)
    if any(abs(visT1_1-i)<0.0001)
        delete(visL1_1);
    end
    if any(abs(visT1_2-i)<0.0001)
        delete(visL1_2);
    end
    if any(abs(visT1_3-i)<0.0001)
        delete(visL1_3);
    end
    if any(abs(visT1_4-i)<0.0001)
        delete(visL1_4);
    end
    if any(abs(visT1_5-i)<0.0001)
        delete(visL1_5);
    end
    if any(abs(visT1_6-i)<0.0001)
        delete(visL1_6);
    end
    if any(abs(visT1_7-i)<0.0001)
        delete(visL1_7);
    end
    if any(abs(visT1_8-i)<0.0001)
        delete(visL1_8);
    end
    if any(abs(visT1_9-i)<0.0001)
        delete(visL1_9);
    end
    if any(abs(visT1_10-i)<0.0001)
        delete(visL1_10);
    end
    if any(abs(visT1_11-i)<0.0001)
        delete(visL1_11);
    end
    if any(abs(visT1_12-i)<0.0001)
        delete(visL1_12);
    end
    if any(abs(visT2_1-i)<0.0001)
        delete(visL2_1);
    end
    if any(abs(visT2_2-i)<0.0001)
        delete(visL2_2);
    end
    if any(abs(visT2_3-i)<0.0001)
        delete(visL2_3);
    end
    if any(abs(visT2_4-i)<0.0001)
        delete(visL2_4);
    end
    if any(abs(visT2_5-i)<0.0001)
        delete(visL2_5);
    end
    if any(abs(visT2_6-i)<0.0001)
        delete(visL2_6);
    end
    if any(abs(visT2_7-i)<0.0001)
        delete(visL2_7);
    end
    if any(abs(visT2_8-i)<0.0001)
        delete(visL2_8);
    end
    if any(abs(visT2_9-i)<0.0001)
        delete(visL2_9);
    end
    if any(abs(visT2_10-i)<0.0001)
        delete(visL2_10);
    end
    if any(abs(visT2_11-i)<0.0001)
        delete(visL2_11);
    end
    if any(abs(visT2_12-i)<0.0001)
        delete(visL2_12);
    end
    if any(abs(visT3_1-i)<0.0001)
        delete(visL3_1);
    end
    if any(abs(visT3_2-i)<0.0001)
        delete(visL3_2);
    end
    if any(abs(visT3_3-i)<0.0001)
        delete(visL3_3);
    end
    if any(abs(visT3_4-i)<0.0001)
        delete(visL3_4);
    end
    if any(abs(visT3_5-i)<0.0001)
        delete(visL3_5);
    end
    if any(abs(visT3_6-i)<0.0001)
        delete(visL3_6);
    end
    if any(abs(visT3_7-i)<0.0001)
        delete(visL3_7);
    end
    if any(abs(visT3_8-i)<0.0001)
        delete(visL3_8);
    end
    if any(abs(visT3_9-i)<0.0001)
        delete(visL3_9);
    end
    if any(abs(visT3_10-i)<0.0001)
        delete(visL3_10);
    end
    if any(abs(visT3_11-i)<0.0001)
        delete(visL3_11);
    end
    if any(abs(visT3_12-i)<0.0001)
        delete(visL3_12);
    end
    if any(abs(visT4_1-i)<0.0001)
        delete(visL4_1);
    end
    if any(abs(visT4_2-i)<0.0001)
        delete(visL4_2);
    end
    if any(abs(visT4_3-i)<0.0001)
        delete(visL4_3);
    end
    if any(abs(visT4_4-i)<0.0001)
        delete(visL4_4);
    end
    if any(abs(visT4_5-i)<0.0001)
        delete(visL4_5);
    end
    if any(abs(visT4_6-i)<0.0001)
        delete(visL4_6);
    end
    if any(abs(visT4_7-i)<0.0001)
        delete(visL4_7);
    end
    if any(abs(visT4_8-i)<0.0001)
        delete(visL4_8);
    end
    if any(abs(visT4_9-i)<0.0001)
        delete(visL4_9);
    end
    if any(abs(visT4_10-i)<0.0001)
        delete(visL4_10);
    end
    if any(abs(visT4_11-i)<0.0001)
        delete(visL4_11);
    end
    if any(abs(visT4_12-i)<0.0001)
        delete(visL4_12);
    end
    if any(abs(visT5_1-i)<0.0001)
        delete(visL5_1);
    end
    if any(abs(visT5_2-i)<0.0001)
        delete(visL5_2);
    end
    if any(abs(visT5_3-i)<0.0001)
        delete(visL5_3);
    end
    if any(abs(visT5_4-i)<0.0001)
        delete(visL5_4);
    end
    if any(abs(visT5_5-i)<0.0001)
        delete(visL5_5);
    end
    if any(abs(visT5_6-i)<0.0001)
        delete(visL5_6);
    end
    if any(abs(visT5_7-i)<0.0001)
        delete(visL5_7);
    end
    if any(abs(visT5_8-i)<0.0001)
        delete(visL5_8);
    end
    if any(abs(visT5_9-i)<0.0001)
        delete(visL5_9);
    end
    if any(abs(visT5_10-i)<0.0001)
        delete(visL5_10);
    end
    if any(abs(visT5_11-i)<0.0001)
        delete(visL5_11);
    end
    if any(abs(visT5_12-i)<0.0001)
        delete(visL5_12);
    end
end
end