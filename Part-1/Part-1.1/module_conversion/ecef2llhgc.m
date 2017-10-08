%% function 'ecef2llh_geocentric'
%
% Transforms coordinates in ECEF to Geocentric LLH
% ECEF: Earth Centered Earth Fixed Frame
% Geocentric LLH: Geocentric Latitude, Longitude, Height Frame
%
% Input  :   pos_ecef           = [x, y, z]' in ECEF frame
% Outputs:   pos_llhgc          = [lat, lon, alt]' in geocentric LLH

% Kuan Chun Hwang
% AERO4701, 2016

function pos_llhgc = ecef2llhgc(pos_ecef)

    % Earth's radius
    global r_earth;

    % Extract data
    x   = pos_ecef(1,:);
    y   = pos_ecef(2,:);
    z   = pos_ecef(3,:);
    
    % Calculating R
    R = sqrt( x.^2 + y.^2 + z.^2);
    
    % Calculating latitude
    lat = asin(z./R);
    
    % Calculating longitude
    lon = atan2( y, x);
    
    % Calculating altitude/height
    h = R - r_earth;
    
    % Output
    pos_llhgc = [lat; lon; h];
end