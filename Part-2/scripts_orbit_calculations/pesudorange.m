%% function pesudorange calculations
%
% Author: Kuan-Chun Hwang

function [x_final,matArray,time_set] = pesudorange(times_pes, pseudor, xsv, ysv, zsv)

% Separating the data sets by time of observation
max_time_pes = length(times_pes);
time_diff = diff(times_pes);
time_diff(max_time_pes) = 1;
t_data_change = find(time_diff == 1);

% Initialize count number
dataNumb = 1;
count = 1;

% Go into a loop that knows how many sets of data there are
for dataSet = 1: length(t_data_change)
    
    % Separate the data sets of each time
    t_set = dataNumb:t_data_change(dataSet);
    
    % we need at least 4 observations
    if length(t_set)>=4

        % Assume x = (0,0,0,0)T initially
        x0_old = [0; 0; 0; 0];
        delta_p0_set = [];
        H_set = [];
        time_set{count} = t_set;
        % Set initial error
        error = 100000;
        
        while max(error) > 0.0000001
            
            % Extract x0 values
            x_o_old = x0_old(1,:);
            y_o_old = x0_old(2,:);
            z_o_old = x0_old(3,:);
            cb_o_old = x0_old(4,:);
            
            % data row count
            i=1;
            
            for n = t_set
                
                % Find the deltap value and H array for each row of data
                [delta_p0,H] = findHandDeltaP(xsv(n),ysv(n),zsv(n),x_o_old,y_o_old,z_o_old,cb_o_old,pseudor(n));
                
                % Extract the set of values
                delta_p0_set(i,1) = delta_p0;
                H_set(i,:) = H;
                i=i+1;
                
            end
            
            deltaX = (inv(H_set'*H_set)) * H_set'* delta_p0_set;
            
            x0_new = x0_old + deltaX;
            
            % Calculate Error
            error = abs(deltaX);
            
            % Update x values
            x0_old = x0_new;
            
        end
        
        % Extract final position and H
        x_final(:,count) = x0_new;
        matArray{count} = H_set;
        count = count + 1;
        
    end
    
    % Move on to next data set
    dataNumb = 1 + t_data_change(dataSet);
    
end

end




