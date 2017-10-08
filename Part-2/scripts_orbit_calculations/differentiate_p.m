%% function differentiate pesudorange equation
%
% Author: Kuan Chun Hwang

function H = differentiate_p(xsv, ysv, zsv, x_o, y_o, z_o)

% From forumula, caculate each row of the H matrix
A = sqrt((xsv - x_o)^2 + (ysv - y_o)^2 + (zsv - z_o)^2);
dp_dx = -(xsv - x_o)/A;
dp_dy = -(ysv - y_o)/A;
dp_dz = -(zsv - z_o)/A;
dp_dcb = 1;

H = [dp_dx, dp_dy, dp_dz, dp_dcb];
end