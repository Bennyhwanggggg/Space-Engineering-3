%% function sate rate equation
%
%  Author: Kuan Chun Hwang

function answer = sate_rate_eqn(x_parameters)

global mu_earth;
global mu_J2_r_sq;

% Extract Parameter values
p = x_parameters(1,:);
f = x_parameters(2,:);
g = x_parameters(3,:);
h = x_parameters(4,:);
k = x_parameters(5,:);
L = x_parameters(6,:);

% Find w,s,r
w = 1 + f.*cos(L) + g.*sin(L);
s_sq = 1 + h.^2 + k.^2;
r = p ./ w;


% Sate Rate Equation
A = [ 0, ((2*p)/w)*(sqrt(p/mu_earth)), 0;
    (sqrt(p/mu_earth))*(sin(L)),(sqrt(p/mu_earth))*(((w+1)*cos(L)+f)/w), -(sqrt(p/mu_earth))*(g/w)*(h*sin(L)-k*cos(L));
    -(sqrt(p/mu_earth))*cos(L),(sqrt(p/mu_earth))*(((w+1)*sin(L)+g)/w), (sqrt(p/mu_earth))*(f/w)*(h*sin(L)-k*cos(L));
    0, 0, (sqrt(p/mu_earth))*((s_sq)*cos(L))/(2*w);
    0, 0, (sqrt(p/mu_earth))*((s_sq)*sin(L))/(2*w);
    0, 0, (sqrt(p/mu_earth))*(1/w)*((h*sin(L))-(k*cos(L)))];

b = [0, 0, 0, 0, 0, (sqrt(mu_earth*p))*((w/p)^2)]';

deltaJ2 = [-(3*mu_J2_r_sq)/(2*(r^4))*(1-(12*(h*sin(L)-k*cos(L))^2)/((1+(h^2)+(k^2))^2));
    (-(12*mu_J2_r_sq)/(r^4)*(((h*sin(L)-k*cos(L))*(h*cos(L)+k*sin(L)))/((1+(h^2)+(k^2))^2)));
    (-(6*mu_J2_r_sq)/(r^4)*(((1-(h^2)-(k^2))*(h*sin(L)-k*cos(L)))/((1+(h^2)+(k^2))^2)))];

fx_dot = A*deltaJ2 + b;

answer = fx_dot;

end