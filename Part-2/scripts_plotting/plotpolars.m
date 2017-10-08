%% function do polar plots
%
% Author Kuan Chun Hwang
function plotpolars(polar_pos_true,polar_pos_est,polar_vis_sat,time_set)

% Get figure and axes
fig = figure;
axe_handle = axes('parent', fig);
% Make sure limits are correct.
lim = polar(0, 90, 'parent', axe_handle);
delete(lim)
% Get the handles of interest (using volatile info above).
handles = findall(fig,'parent', axe_handle, 'Type', 'text');
handles = handles(strncmp('  ', get(handles,'String'), 2));
handles = sort(handles);
% Relabel from inside out.
labels = {'72','54', '36', '18', '0'};
for i = 1:5
    
    set(handles(i),'String', labels{i})
    
end
hold on

% Plot polar plots
polar(polar_pos_true(2,:),abs(90-rad2deg(polar_pos_true(3,:))),'r','parent', axe_handle);

hold on
polar(polar_pos_est(2,:),abs(90-rad2deg(polar_pos_est(3,:))),'b','parent', axe_handle)

hold on

% Rotate plots
view([90 -90])
ti = title('UAVs trajectory polar plot w.r.t ground station');
pos = get(ti,'Position');
set(ti,'Position',[pos(1)+8 pos(2) pos(3)]);
set(ti,'FontSize',8);
legend('True','Estimated')
% Animation
for i = 1:length(polar_pos_est)
    est_trace =  polar(polar_pos_est(2,i),abs(90-rad2deg(polar_pos_est(3,i))),'o','parent', axe_handle);
    vis_sat_plot = polar(polar_vis_sat(2,time_set{i}),abs(90-rad2deg(polar_vis_sat(3,time_set{i}))),'.','parent', axe_handle);
    set(vis_sat_plot,'markersize',12)
    pause (0.001);
    if i~= length(polar_pos_est)
        delete(est_trace);
        delete(vis_sat_plot);
    end
    
end

end