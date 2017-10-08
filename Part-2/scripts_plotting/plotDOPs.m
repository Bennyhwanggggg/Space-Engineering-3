%% function plot DOPS and also determine the time when best and worst configuration occur
%
% Author Kuan Chun Hwang

function [maxDOP,minDOP]=plotDOPs(DOPs,rep_t)

% Extract values
GDOP = DOPs(1,:);
PDOP = DOPs(2,:);
HDOP = DOPs(3,:);
VDOP = DOPs(4,:);
TDOP = DOPs(5,:);

% Select a DOP to plot
prompt = {'Choose GDOP, PDOP, HDOP, VDOP, TDOP'};
dlg_title = 'Select DOP to plot';
answer = inputdlg(prompt,dlg_title);

tr = strcmp('GDOP',answer);
if tr == 1
    figure
    plot(rep_t,GDOP)
    xlabel('time(s)')
    ylabel('GDOP')
    title('GDOP')
    
    % Find the Higest and Lowest GDOP
    maxDOP = find(GDOP == max(GDOP));
    minDOP = find(GDOP == min(GDOP));
end

tr = strcmp('PDOP',answer);
if tr == 1
    figure
    plot(rep_t,PDOP)
    xlabel('time(s)')
    ylabel('PDOP')
    title('PDOP')
    
    % Find the Higest and Lowest PDOP
    maxDOP = find(PDOP == max(PDOP));
    minDOP = find(PDOP == min(PDOP));
end

tr = strcmp('HDOP',answer);
if tr == 1
    figure
    plot(rep_t,HDOP)
    xlabel('time(s)')
    ylabel('HDOP')
    title('HDOP')
    
    % Find the Higest and Lowest HDOP
    maxDOP = find(HDOP == max(HDOP));
    minDOP = find(HDOP == min(HDOP));
    
end

tr = strcmp('VDOP',answer);
if tr == 1
    figure
    plot(rep_t,VDOP)
    xlabel('time(s)')
    ylabel('VDOP')
    title('VDOP')
    
    % Find the Higest and Lowest VDOP
    maxDOP = find(VDOP == max(VDOP));
    minDOP = find(VDOP == min(VDOP));
end

tr = strcmp('TDOP',answer);
if tr == 1
    figure
    plot(rep_t,TDOP)
    xlabel('time(s)')
    ylabel('TDOP')
    title('TDOP')
    
    % Find the Higest and Lowest TDOP
    maxDOP = find(TDOP == max(TDOP));
    minDOP = find(TDOP == min(TDOP));
end

end
