%% function calculate period of an orbit
%
% Author Kuan Chun Hwang

function period = calc_period(a)

global mu_earth

MeanMotion  = sqrt(mu_earth/a^3); % Mean Motion (rad/s)
period = 2*pi/MeanMotion;  % Period (s)

end