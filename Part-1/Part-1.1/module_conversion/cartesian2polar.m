%% function 'cartesian2polar'
%
% Transforms coordinates in cartesian coordinates [x, y, z]' 
% into polar coordinates [R, az, el]'
%
% Input  : pos_cartesian = [x,  y,  z]' in cartesian coordinates [m, m, m]
% Output : pos_polar     = [R, az, el]' in polar coordinates [m, rad, rad]
%
% Kuan Chun Hwang


function pos_polar = cartesian2polar(pos_cartesian)

% Extract Data
x = pos_cartesian(1,:);
y = pos_cartesian(2,:);
z = pos_cartesian(3,:);

% Calculate Radius
R=sqrt( x.^2 + y.^2 + z.^2);

% Calculate az
az = atan( y ./ x);

% Calculate el
el = atan2( -z , sqrt(x.^2+y.^2));

[ pos_polar ]  = [R; az; el];    

end