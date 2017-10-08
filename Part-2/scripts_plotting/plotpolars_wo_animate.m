%% function do polar plots without animation
%
% Author Kuan Chun Hwang
function plotpolars_wo_animate(polar_pos_true,polar_pos_est)

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
ti = title('Satellite polar plot w.r.t Ground Station');
pos = get(ti,'Position');
set(ti,'Position',[pos(1)+8 pos(2) pos(3)]);
set(ti,'FontSize',8);
legend('initial','final')
