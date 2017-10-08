%% function find ecef position of observation satellites
%
% Author: Kuan Chun Hwang

function [xsv, ysv, zsv] = find_ecef_obs(n_sat,a,e,i_deg, cap_omega_deg, omega_deg, M_deg,Epoc_time, times_pes, n_sat_pes,t_vernal)

[pos_eci_ephem, pos_llhgc_ephem, pos_ecef] = keplerian(n_sat,a,e,i_deg, cap_omega_deg, omega_deg, M_deg,Epoc_time,times_pes',t_vernal);
% Extract pseudorange data
% Satellite number in pesudorange data

for n = 1:length(n_sat_pes)
    
    % Extract satellite number
    sat_n_pes = n_sat_pes(n);

    % Extract satellite ecef data
    xsv(n) = pos_ecef(n,sat_n_pes,1);
    ysv(n) = pos_ecef(n,sat_n_pes,2);
    zsv(n) = pos_ecef(n,sat_n_pes,3);
end

end