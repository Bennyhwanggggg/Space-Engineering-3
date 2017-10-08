%% function "Plot ground station observation range"
%
%  Author: Kuan- Chun Hwang

function plot_observations(pos_llhgc,visible_time,gs_number,n_sat,gs_llh)

% Convert from radians to degree
gs_lat_deg = rad2deg(gs_llh(1,:));
gs_lon_deg = rad2deg(gs_llh(2,:));

% Ground plot
figure
map = imread('Map.jpg');
image([-180 180], [90 -90], map);
axis xy
hold on
plot(pos_llhgc(:, :, 2), pos_llhgc(:, :, 1));

hold on

for gs_n = 1:length(gs_number)
    for n = 1:length(n_sat)
        
        % Plot the observation area on graph
        scatter(pos_llhgc(visible_time{gs_n,n}, n, 2), pos_llhgc(visible_time{gs_n,n}, n, 1),'r','LineWidth',0.1);
        hold on
    end
end

% Plot Ground Station Locations
scatter(gs_lon_deg,gs_lat_deg,'g','LineWIdth',2)

title('Ground Station Observation Range');
xlabel('longitude');
ylabel('latitude');
end