%% function 'llh_geocentric2ecef'
%
% Transforms coordinates in Geocentric LLH to coordinates in ECEF
% Geocentric LLH: Geocentric Latitude, Longitude, Height Frame
% ECEF: Earth Centered Earth Fixed Frame
%
% Input  : pos_llhgc          = [lat, lon, alt]' in geocentric LLH
% Output : pos_ecef           = [x, y, z]' in ECEF frame
%
% Kuan Chun Hwang

function pos_ecef = llhgc2ecef(pos_llhgc)
    
    % Earth's radius
    global r_earth;

    % Extract Data
    lat = pos_llhgc(1,:);
    lon = pos_llhgc(2,:);
    h = pos_llhgc(3,:);
    
    % Output
    pos_ecef = [ (r_earth+h).*cos(lat).*cos(lon); (r_earth+h).*cos(lat).*sin(lon); (r_earth+h).*sin(lat)];
    
end