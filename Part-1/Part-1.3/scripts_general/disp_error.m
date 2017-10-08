%% function display error data
%
% Author: Kuan Chun Hwang

function disp_error(differnce)

% Display the results
fprintf('Semi-Major axis error in percentage at 200db = %g\n', differnce(1));
fprintf('Inclination error in percentage at 200db = %g\n', differnce(2));
fprintf('Eccentrcity error at 200db = %g\n', differnce(3));
fprintf('Right Ascension of Ascending Node error in percentage at 200db = %g\n', differnce(4));
fprintf('Argument of Perigee error in percentage at 200db = %g\n', differnce(5));
fprintf('Mean anomaly error in percentage at 200db = %g\n', differnce(6));

end