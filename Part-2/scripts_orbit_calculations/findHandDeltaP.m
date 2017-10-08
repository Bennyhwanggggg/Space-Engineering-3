%% function find deltaP0 and H
%
% Auhtor Kuan Chun Hwang

% Using the ECEF coordinates of the same n from previously
% calculations
function [delta_p0,H] = findHandDeltaP(xsv,ysv,zsv,x_o_old,y_o_old,z_o_old,cb_o_old,pseudor)

% Using the ECEF coordinates of the same n from previously
% calculations
fx = sqrt(((xsv-x_o_old)^2) + ((ysv - y_o_old)^2) + ((zsv - z_o_old)^2)) + cb_o_old;

delta_p0 = pseudor - fx;

% H = d(f(x)/dx) a tx= x0
H = differentiate_p(xsv, ysv, zsv, x_o_old, y_o_old, z_o_old);
end