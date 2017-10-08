%% Calculate deltaV in ECI frame
% This function converts velocity vectors from local frame to ECI
%
% Author Kuan Chun Hwang

function deltaV_ECI   = deltaV(deltaV_normalCoord,v,r,theta,psi)

% Compute x direction unit vector
x_local = v/norm(v);

% Calculate the y direction vector
y_local = cross(v,r);
% Compute y direction unit vector
y_local = y_local/norm(y_local);

% Calculate the y direction vector
z_local = cross(x_local,y_local);
% Compute z direction unit vector
z_local = z_local/norm(z_local);

% Compute orthogonal matrix that defines principle axes of an inertial velocity coordinate frame.
local2eci = [x_local,y_local,z_local];

% Find change in velocity in Local frame
deltaV = deltaV_normalCoord*[cos(theta)*cos(psi); cos(theta)*sin(psi); sin(theta)];

% Convert change in velocity from local frame to ECI frame
deltaV_ECI = local2eci*deltaV;

end