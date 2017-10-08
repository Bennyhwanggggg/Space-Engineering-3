%% function obtain satellite sensor readings for each ground station
% Author Kuan Chun Hwang

function [sat_lg_polar_recieved, satChosen_vis_t, sat_lg_polar] = obtainReadings(satChosen,visible_time,gs_number,gs_llh,sat_pos_ecef)

% Find satellite ecef
gs_ecef = llhgc2ecef(gs_llh);

% Retrieve the satellite's LG position w.r.t to each ground station during
% their respective visible time
for gs_n = gs_number
    
    for n = satChosen
        
        % Retrieve the visible time period of the chosen satellite for a ground
        % station
        satChosen_vis_t = visible_time{gs_n,n};
        
        % Extract satellite ecef position
        sat_ecef = [(sat_pos_ecef(satChosen_vis_t,n,1))';(sat_pos_ecef(satChosen_vis_t,n,2))';(sat_pos_ecef(satChosen_vis_t,n,3))'];
        
        % Caculate satellites relative position to ground station in ecef
        sat_ecef_wrt_gs = bsxfun(@minus,sat_ecef,gs_ecef(:,gs_n));
        
        % Convert satellite relative position to GS from ECEF to LGCV
        sat_lg_cart = ecef2lg(sat_ecef_wrt_gs,gs_llh(:,gs_n));
        
        % Convert from cartesian to polar
        sat_lg_polar = cartesian2polar(sat_lg_cart);
        
        % For a simple model, we will add random error to the simulated
        % truth. 
        range_error = -100000+(100000+100000).*rand(1,length(sat_lg_polar));
        az_el_error = deg2rad(-20)+(deg2rad(20)-deg2rad(-20)).*rand(2,length(sat_lg_polar));        
        signal_error = [range_error;az_el_error];
        
        % Add noise error on to readings
        sat_lg_polar_recieved{gs_n,n} = sat_lg_polar + signal_error;
        
    end
end

end