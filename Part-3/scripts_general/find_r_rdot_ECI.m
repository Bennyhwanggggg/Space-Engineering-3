%% Calculate r and r_dot ECI vectors
%
% Author Kuan Chun Hwang

function [r, rDot] = find_r_rdot_ECI(pos_eci,e,i_rad,omega_rad,RAAN_rad,true_anomaly_rad,p)

% Initialize global constants
global mu_earth

% Extract ECI coordinates
rX = pos_eci(1,1,1);
rY = pos_eci(1,1,2);
rZ = pos_eci(1,1,3);

r = [rX; rY; rZ];

% Find velocity in ECI
rDotX_orbital = -sqrt(mu_earth/p)*sin(true_anomaly_rad);
rDotY_orbital =  sqrt(mu_earth/p)*(e+cos(true_anomaly_rad));
rDotZ_orbital  = 0;

orbit_to_ECI = orbit2ECI(RAAN_rad,omega_rad,i_rad);

rDot = orbit_to_ECI*[rDotX_orbital;rDotY_orbital;rDotZ_orbital];

end