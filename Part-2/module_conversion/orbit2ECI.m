%% function convert from perifocal to ECI
%
% Author Kuan Chun Hwang

function orbit_to_ECI = orbit2ECI(cap_omega_radn,omega_radn,i_radn)

% Initialize Transfer matrix from orbital coordinates to ECI frame
orbit_to_ECI = [cos(cap_omega_radn)*cos(omega_radn)-sin(cap_omega_radn)*sin(omega_radn)*cos(i_radn) -cos(cap_omega_radn)*sin(omega_radn)-sin(cap_omega_radn)*cos(omega_radn)*cos(i_radn) sin(cap_omega_radn)*sin(i_radn)
    sin(cap_omega_radn)*cos(omega_radn)+cos(cap_omega_radn)*sin(omega_radn)*cos(i_radn) -sin(cap_omega_radn)*sin(omega_radn)+cos(cap_omega_radn)*cos(omega_radn)*cos(i_radn) -cos(cap_omega_radn)*sin(i_radn)
    sin(omega_radn)*sin(i_radn) cos(omega_radn)*sin(i_radn) cos(i_radn)];

end