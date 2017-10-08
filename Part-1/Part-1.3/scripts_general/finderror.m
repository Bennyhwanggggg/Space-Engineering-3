%% function find error
%
% Author: Kuan Chun Hwang

function noise_effect = finderror(orig_para,orb_parameters)

% Initialize global variable
global mu_earth;

% Extract Data of TLE
inclination_i = orig_para(1);
e_i = orig_para(2);
capital_omega_i = orig_para(3);
omega_i = orig_para(4);
M_i = orig_para(5);

% Find initial semi major axis
a_i = nthroot(mu_earth/(orig_para(6)^2),3); %(m)

original = [a_i;inclination_i;e_i;capital_omega_i;omega_i;M_i];

% Extract Data from tracking 
a = orb_parameters(1);
inclination = rad2deg(orb_parameters(2));
e = orb_parameters(3);
capital_omega = rad2deg(orb_parameters(4));
omega = rad2deg(orb_parameters(5));
theta = orb_parameters(6);

% Find mean anomaly
essentric_anomaly = 2*atan((tan(theta/2))/(sqrt((1+e)/(1-e))));
M = rad2deg((essentric_anomaly - e*sin(essentric_anomaly)));

Tracking = [a;inclination;e;capital_omega;omega;M];

% Find difference
diff = Tracking - original;

noise_effect = [abs(100*diff(1)/original(1)); abs(100*diff(2)/original(2)); diff(3); abs(100*diff(4)/original(4)); abs(100*diff(5)/original(5)); abs(100*diff(6)/original(6))];

end


