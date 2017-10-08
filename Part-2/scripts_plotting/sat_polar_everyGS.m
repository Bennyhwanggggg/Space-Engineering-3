%% function "Calculate number of observations from each ground station"
%
%  Author: Kuan- Chun Hwang

function sat_polar_everyGS(gs_number,gs_llh,satChosen,sat_pos_ecef_true,sat_pos_ecef_est)

% Find satellite ecef
gs_ecef = llhgc2ecef(gs_llh);

for gs_n = 1:length(gs_number)
    
    for n = 1:length(satChosen)
        
        % Extract satellite ecef position
        sat_ecef = [(sat_pos_ecef_true(:,n,1))';(sat_pos_ecef_true(:,n,2))';(sat_pos_ecef_true(:,n,3))'];
        
        % Caculate satellites relative position to ground station in ecef
        sat_ecef_wrt_gs = bsxfun(@minus,sat_ecef,gs_ecef(:,gs_n));
        
        % Convert satellite relative position to GS from ECEF to LGCV
        sat_lg_cart = ecef2lg(sat_ecef_wrt_gs,gs_llh(:,gs_n));
        
        % Convert from cartesian to polar
        sat_lg_polar_i = cartesian2polar(sat_lg_cart);
        
        % Extract satellite ecef position
        sat_ecef = [(sat_pos_ecef_est(:,n,1))';(sat_pos_ecef_est(:,n,2))';(sat_pos_ecef_est(:,n,3))'];
        
        % Caculate satellites relative position to ground station in ecef
        sat_ecef_wrt_gs = bsxfun(@minus,sat_ecef,gs_ecef(:,gs_n));
        
        % Convert satellite relative position to GS from ECEF to LGCV
        sat_lg_cart = ecef2lg(sat_ecef_wrt_gs,gs_llh(:,gs_n));
        
        % Convert from cartesian to polar
        sat_lg_polar_f = cartesian2polar(sat_lg_cart);
        
        plotpolars_wo_animate(sat_lg_polar_i,sat_lg_polar_f)
        
    end
    
end

end
