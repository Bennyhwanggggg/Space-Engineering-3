%% Store pseudorange data into data sets and put it into time, number of satellites and pseudorange.
%
% Author Kuan Chun Hwang

function [times_pse, pseudor_pse, xsv, ysv, zsv]= error_sat_findFix(times_pse,n_sat_pse, pseudor_pse, err_index, xsv, ysv, zsv)

% Separating the data sets by time of observation
max_time_pes = length(times_pse);
time_diff = diff(times_pse);
time_diff(max_time_pes) = 1;
t_data_change = find(time_diff == 1);

% Initialize count number
dataNumb = 1;
count = 1;

for dataSet = 1: length(t_data_change)
    
    if dataSet > min(err_index) && dataSet < max(err_index)
        
        % Separate the data sets of each time
        for t_set = dataNumb:t_data_change(dataSet);
            
            sat_n_use(count,1) = n_sat_pse(t_set);
            count = count +1;
        
        end
        
    end
    
end

% Find the satellite that appears the most
sat_n_most = mode(sat_n_use);

for dataSet = 1: length(t_data_change)  
    
    if dataSet > min(err_index) && dataSet < max(err_index)
        
        % Separate the data sets of each time
        for t_set = dataNumb:t_data_change(dataSet);
            
            if n_sat_pse(t_set) == sat_n_most
                
                n_sat_pse(t_set) = [];
                times_pse(t_set) = [];
                pseudor_pse(t_set) = [];
                xsv(t_set) = [];
                ysv(t_set) = [];
                zsv(t_set) = [];
            
            end
            
        end
        
    end
end
end