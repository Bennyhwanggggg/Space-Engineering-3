%% function 'llh_geocentric2ecef'
%
% Transforms coordinates in Geodetic LLH to coordinates in ECEF
% Geodetic LLH: Geodetic Latitude, Longitude, Height Frame
% ECEF: Earth Centered Earth Fixed Frame
%
% Input  : pos_llhgd          = [lat, lon, alt]' in geodetic LLH
% Output : pos_ecef           = [x, y, z]' in ECEF frame
%
% Kuan Chun Hwang

function pos_ecef = llhgd2ecef(pos_llhgd)
   
    % Semi-Major Axis of Earth
    global r_earth;
    
    % Eccentricity of Earth
    global e_earth;

    % Extract Data
    lat = pos_llhgd(1,:);
    lon = pos_llhgd(2,:);
    h = pos_llhgd(3,:);
    
    % Finding N
    N = r_earth./sqrt(1 - e_earth^2.*(sin(lat).^2));
    
    % Output
    pos_ecef = [(N+h).*cos(lat).*cos(lon); (N+h).*cos(lat).*sin(lon); ((N.*(1-e_earth^2)+h).*sin(lat))];
    
end