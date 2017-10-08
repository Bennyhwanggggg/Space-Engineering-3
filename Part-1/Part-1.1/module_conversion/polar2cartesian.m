%% function 'polar2cartesian'
% Transforms polar coordinates [R, az, el]' into
% cartesian coordinates [x, y, z]' 

% Input  : pos_polar     = [R, az, el]' in polar coordinates [m, rad, rad]
% Output : pos_cartesian = [x,  y,  z]' in cartesian coordinates [m, m, m]
%
% Kuan Chun Hwang


function pos_cartesian = polar2cartesian(pos_polar)

% Extract Data
R = pos_polar(1,:);
az = pos_polar(2,:);
el = pos_polar(3,:);

% Output
pos_cartesian = [ R*cos(el).*cos(az); R.*cos(el).*sin(az); -R.*sin(el)];

end