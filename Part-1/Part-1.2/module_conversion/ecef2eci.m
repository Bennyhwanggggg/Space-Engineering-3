%% function 'ecef2eci'
%
% Transforms coordinates in ECEF to coordinates in ECI
% ECEF: Earth Centered Earth Fixed Frame
% ECI : Earth Centered Inertial Frame
%
% Input  : pos_ecef = [x; y; z] | ECEF                     [m]
%          times    = times since vernal equinox alignment [s]
% Output : pos_eci  = [x; y; z] | ECI                      [m]
%
% Kuan Chun Hwang


function pos_eci = ecef2eci(pos_ecef, times)

    % This is the rotation rate of Earth (rad/s)
    global w_earth;
    
    % Transformation matrix from ECI to ECEF
    [ ECI_to_ECEF ] = [ cos(w_earth*times) sin(w_earth*times) 0; -sin(w_earth*times) cos(w_earth*times) 0; 0 0 1];

    % Transformation matrix from ECEF to ECI
    [ ECEF_to_ECI ] = transpose(ECI_to_ECEF);

    % Output
    [ pos_eci ] = ECEF_to_ECI * pos_ecef;

end