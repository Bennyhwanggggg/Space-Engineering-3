%% function 'lg2ecef'
%
% Transforms coordinates in LG to coordinates in ECEF
% LG  : Local Geocentric Vertical Frame OR Local Geodetic Vertical Frame
% ECEF: Earth Centered Earth Fixed Frame
%
% Input  : pos_lg                    = [x, y, z]' in LG frame
%          pos_llh_ground            = [lat, lon, h]'|Observer
% Output : pos_ecef                  = [x, y, z]' in ECEF frame
%
% Kuan Chun Hwang

function pos_ecef = lg2ecef(pos_lg, pos_llh_ground)

% Extract Data
lat = pos_llh_ground(1,:);
lon = pos_llh_ground(2,:);

% Transformation matrix from LGDV to ECEF
[ LGDV_to_ECEF ] = [ -sin(lat).*cos(lon) -sin(lon) -cos(lat).*cos(lon); -sin(lat).*sin(lon) cos(lon) -cos(lat).*sin(lon); cos(lat) 0 -sin(lat)];

% Output
[ pos_ecef ] = LGDV_to_ECEF * pos_lg;


end