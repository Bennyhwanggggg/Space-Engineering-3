%% function DOP best and worse config plot
%
% Author Kuan Chun Hwang

function DOPconfig_plot(maxDOP,minDOP,polar_vis_sat,time_set)

% Get figure and axes
fig = figure;
axe_hand = axes('parent', fig);
% Make sure limits are correct.
lim = polar(0, 90, 'parent', axe_hand);
delete(lim)
% Get the handles of interest (using volatile info above).
handles = findall(fig,'parent', axe_hand, 'Type', 'text');
handles = handles(strncmp('  ', get(handles,'String'), 2));
handles = sort(handles);
% Relabel from inside out.
labels = {'72','54', '36', '18', '0'};
for i = 1:5
    
    set(handles(i),'String', labels{i})
    
end
hold on

polar(polar_vis_sat(2,time_set{maxDOP}),abs(90-rad2deg(polar_vis_sat(3,time_set{maxDOP}))),'.','parent', axe_hand);

% Rotate plots
view([90 -90])
ti = title('Worst Satellite Configuration');
pos = get(ti,'Position');
set(ti,'Position',[pos(1)+8 pos(2) pos(3)]);
set(ti,'FontSize',8);

% Get figure and axes
fig = figure;
axe_hand = axes('parent', fig);
% Make sure limits are correct.
lim = polar(0, 90, 'parent', axe_hand);
delete(lim)
% Get the handles of interest (using volatile info above).
handles = findall(fig,'parent', axe_hand, 'Type', 'text');
handles = handles(strncmp('  ', get(handles,'String'), 2));
handles = sort(handles);
% Relabel from inside out.
labels = {'72','54', '36', '18', '0'};
for i = 1:5
    
    set(handles(i),'String', labels{i})
    
end
hold on

polar(polar_vis_sat(2,time_set{minDOP}),abs(90-rad2deg(polar_vis_sat(3,time_set{minDOP}))),'.','parent', axe_hand);

% Rotate plots
view([90 -90])
ti = title('Best Satellite Configuration');
pos = get(ti,'Position');
set(ti,'Position',[pos(1)+8 pos(2) pos(3)]);
set(ti,'FontSize',8);
