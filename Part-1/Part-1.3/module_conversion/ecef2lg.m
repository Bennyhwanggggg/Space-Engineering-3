%% function 'ecef2lg'
%
% Transforms coordinates in ECEF to coordinates in LG
% ECEF: Earth Centered Earth Fixed Frame
% LG  : Local Geocentric Vertical Frame OR Local Geodetic Vertical Frame
%
% Input  : pos_ecef                  = [x, y, z]' in ECEF frame
%          pos_llh_ground            = [lat, lon, h]'|Observer
%
% Outputs: pos_lg                    = [x, y, z]' in LG frame
%
% Kuan Chun Hwang


function pos_lg = ecef2lg(pos_ecef, pos_llh_ground)
    
% Extract Data
lat = pos_llh_ground(1,:);
lon = pos_llh_ground(2,:);

% Transformation matrix from LGDV to ECEF
[ LG_to_ECEF ] = [ -sin(lat).*cos(lon) -sin(lon) -cos(lat).*cos(lon); -sin(lat).*sin(lon) cos(lon) -cos(lat).*sin(lon); cos(lat) 0 -sin(lat)];

% Transformation matrix from ECEF to LGDV
[ ECEF_to_LG ] = transpose(LG_to_ECEF);

% Output
[ pos_lg ]= ECEF_to_LG * pos_ecef;

end