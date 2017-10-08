%% function find used times in pesudorange data
%
% Author Kuan Chun Hwang

function [rep_t, vis_sat_n,t_obs] = find_t_repeat(obs_time)

% Seprate the time into data sets
max_time_pes = length(obs_time);
time_diff = diff(obs_time);
time_diff(max_time_pes) = 1;
t_data_change = find(time_diff == 1);

% Initialize count numbers
dataNumb = 1;
count = 1;

for dataSet = 1: length(t_data_change)
    
    % Separate the data sets of each time
    t_set = dataNumb:t_data_change(dataSet);
    vis_sat_n(dataSet) = length(t_set);
    t_obs(dataSet) = obs_time(dataNumb);
    
    % we need at least 4 observations for it to be useful
    if length(t_set)>=4
        
        % Extract the time value
        rep_t(count) = obs_time(dataNumb);
        count = count + 1;
    end

    % Move on to next data set
    dataNumb = 1 + t_data_change(dataSet);
end

end