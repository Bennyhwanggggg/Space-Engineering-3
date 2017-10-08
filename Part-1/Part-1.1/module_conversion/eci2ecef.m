%% function 'eci2ecef_multitime'
%
% Transforms coordinates in ECI to coordinates in ECEF
% ECI : Earth Centered Inertial Frame
% ECEF: Earth Centered Earth Fixed Frame
%
% Input  : pos_eci  = [x; y; z] | ECI                      [m]
%          times    = times since vernal equinox alignment [s]
% Output : pos_ecef = [x; y; z] | ECEF                     [m]
%
% Kuan Chun Hwang


function pos_ecef = eci2ecef(pos_eci, times)

    % This is the rotation rate of Earth (rad/s)
    global w_earth;
    
    % Transformation matrix from ECI to ECEF
    [ ECI_to_ECEF ] = [ cos(w_earth*times) sin(w_earth*times) 0; -sin(w_earth*times) cos(w_earth*times) 0; 0 0 1];
    
    % Output   
    [ pos_ecef ] =  ECI_to_ECEF * pos_eci;

end