%% function find observation jacobain matrix dy/dx
%
% Author Kuan Chun Hwang

function dy_dx = ObservationJacobian(est_lg_cart, gs_llh_use, dt)

% Initiate global variables
global w_earth;

% Extract xLGV yLGV zLGV and dt
x = est_lg_cart(1,:);
y = est_lg_cart(2,:);
z = est_lg_cart(3,:);

% Times is the time since vernal equinox
times = dt;

xyzSQ = (x^2)+(y^2)+(z^2);
xySQ = (x^2)+(y^2);

% Solve for relevant variables
dp_dx = x / sqrt(xyzSQ);
dp_dy = y / sqrt(xyzSQ);
dp_dz = z / sqrt(xyzSQ);

daz_dx = -y / xySQ;
daz_dy = x / xySQ;
daz_dz = 0;

dtheta_dx = (x*z) / (sqrt(xySQ)*xyzSQ);
dtheta_dy = (y*z) / (sqrt(xySQ)*xyzSQ);
dtheta_dz = (-sqrt(xySQ)/xyzSQ);

% Compute dy/drLGV
dy_drLGV = [dp_dx dp_dy dp_dz; daz_dx daz_dy daz_dz; dtheta_dx dtheta_dy dtheta_dz];

% Extract latitude and longitude
lat = gs_llh_use(1,:);
lon = gs_llh_use(2,:);

% Transformation matrix from LGDV to ECEF
[ LG_to_ECEF ] = [ -sin(lat).*cos(lon) -sin(lon) -cos(lat).*cos(lon); -sin(lat).*sin(lon) cos(lon) -cos(lat).*sin(lon); cos(lat) 0 -sin(lat)];

% Transformation matrix from ECEF to LGDV
[ ECEF_to_LG ] = transpose(LG_to_ECEF);

% Transformation matrix from ECI to ECEF
[ ECI_to_ECEF ] = [ cos(w_earth*times) sin(w_earth*times) 0; -sin(w_earth*times) cos(w_earth*times) 0; 0 0 1];

% Compute drLGV/dx
drLGV_dx = [ECEF_to_LG*ECI_to_ECEF zeros(3)];

% Compute dy/dx
dy_dx = dy_drLGV*drLGV_dx;

end

