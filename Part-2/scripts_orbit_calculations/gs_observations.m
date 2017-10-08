%% function "Calculate number of observations from each ground station"
%
%  Author: Kuan- Chun Hwang

function [percentage_obs,visible_time] =  gs_observations(gs_number,gs_llh,n_sat,sat_pos_ecef,times)

% Find satellite ecef
gs_ecef = llhgc2ecef(gs_llh);

% Intialize visible time matrix
total_gs_n = length(gs_number);
total_sat_n = length(n_sat);
visible_time = {total_gs_n,total_sat_n};
% perct_t_visible = zeros(total_gs_n,total_sat_n,1);

count = 1;
for gs_n = 1:length(gs_number)
    for n = 1:length(n_sat)
            
            % Extract satellite ecef position
            sat_ecef = [(sat_pos_ecef(:,n,1))';(sat_pos_ecef(:,n,2))';(sat_pos_ecef(:,n,3))'];
            
            % Caculate satellites relative position to ground station in ecef
            sat_ecef_wrt_gs = bsxfun(@minus,sat_ecef,gs_ecef(:,gs_n));
            
            % Convert satellite relative position to GS from ECEF to LGCV
            sat_lg_cart = ecef2lg(sat_ecef_wrt_gs,gs_llh(:,gs_n));
            
            % Convert from cartesian to polar
            sat_lg_polar = cartesian2polar(sat_lg_cart);
            
            % Store visible time into cells
            visible_t{count} = find(sat_lg_polar(3,:)>0);
            visible_time{gs_n,n} = visible_t{count};
            
            % find percentage of time visible for the specfic satellite
            % from that ground station
            percentage_obs{gs_n,n} = (length(visible_time{gs_n,n})/length(times))*100;
            
            % Update count
            count = count + 1;
    end
    
end

% Store percentage of visible time into matrix
perct_t_visible(:,:,:) = reshape(percentage_obs,total_gs_n,total_sat_n ,1);

f = figure('Name','time% a satellite in been observed by a ground station (GS_n vs Sat_n)','NumberTitle','off');
 
% Create the uitable
t = uitable(f,'Data',perct_t_visible);

end
