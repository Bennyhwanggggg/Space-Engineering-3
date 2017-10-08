%% function Universal Conic soultion
%
% Author Kuan Chun Hwang

function [x0] = UCM(r2,v2,time_use2,gs_use,gs_llh,sat_lg_polar_reci_noise,times_obs)

% Initialize global variable
global mu_earth;

% Get ground station  llh
gs_llh_use = gs_llh(:,gs_use);
gs_ecef = llhgc2ecef(gs_llh_use);

% Set initial parameters to be estimated
r_t0 = r2;
r_dot_t0 = v2;
R_t0 = norm(r_t0);
R_dot_t0 = norm(r_dot_t0);

% Initial time. This initial time is the fourth timestep of the visible
% time period
t0 = time_use2;

times_obs(1:4) = [];
% Set initial error
deltaX = 1000;
iter = 1;

% Start iteration
while max(abs(deltaX)) > 0.001
    
    % Compute semi-Major axis using vis viva law
    a = ((2/R_t0) - ((R_dot_t0^2)/mu_earth))^-1;
    
    % Set initial phi and other inital variables
    n = 1;
    sat_t = 5;
    phi = pi;
    
    for t = times_obs
        
        % Find t - t0
        dt = t - t0;
        
        % Compute change in eccentric anomaly
        phi_error = 1000;
        while phi_error > 0.001
            fE = sqrt((a^3)/mu_earth)*(phi - (1 - (R_t0/a))*sin(phi) + (r_t0'*r_dot_t0)/sqrt(mu_earth*a)*(1-cos(phi))) - dt;
            dE = sqrt((a^3)/mu_earth)*(1 - (1 - (R_t0/a))*cos(phi) + (r_t0'*r_dot_t0)/sqrt(mu_earth*a)*sin(phi));
            
            phiNew = phi - fE/dE;
            
            phi_error = abs(phiNew - phi);
            
            phi = phiNew;
        end
        
        % Find r(t) and velocity
        [f,g,f_dot,g_dot,r,r_t,r_dot_t] = find_r_t_and_r_dot_t(a,phi,dt,r_t0,r_dot_t0,R_t0,R_dot_t0);
        
        % Compute magnitude of r_t and r_dot_t
        R_t = norm(r_t);
        R_dot_t = norm(r_dot_t);
        
        % Compute state transition Jacobian Matrix
        dxt_dxt0 = StateTransition(dt, a, f, g, f_dot, g_dot, r_t0, r_dot_t0, R_t0, R_dot_t0, r_t, r_dot_t, R_t, R_dot_t);
        
        % Compute observation model Jacobian
        % Find estimated position in LGCV first in polar coordinates and
        % cartesian coordinates
        est_eci = r_t;
        est_ecef = eci2ecef(est_eci,t);
        est_ecef_wrt_gs = est_ecef - gs_ecef;
        est_lg_cart = ecef2lg(est_ecef_wrt_gs,gs_llh_use);
        est_lg_pol = cartesian2polar(est_lg_cart);
        
        yhat = est_lg_pol;
        
        dy_dx = ObservationJacobian(est_lg_cart, gs_llh_use, t);
        
        % Compute H, Jacobian Matrix
        H(n:n+2,:) = dy_dx * dxt_dxt0;
        
        % Import true position
        y_true = sat_lg_polar_reci_noise(:,sat_t);
        
        % Compute DeltaY
        deltaY(n:n+2,:) = y_true - yhat;
        
        % Update index variables
        n = n+3;
        sat_t = sat_t + 1;
        
    end
    
    % Compute deltaX
    deltaX = ((H'*H)\ H') * deltaY;
    
    % Store deltaX into another variable for plotting purpose
    dX(iter) = max(deltaX);
    
    r_t0 = r_t0 + deltaX(1:3);
    r_dot_t0 = r_dot_t0 + deltaX(4:6);
    R_t0 = norm(r_t0);
    R_dot_t0 = norm(r_dot_t0);
    
    % Update index variable
    iter = iter+1;
end

% Extract converged result
x0 = [r_t0 r_dot_t0 ];

% Shoe plot of convergence
figure
plot(1:(iter-1),dX);
title('convergence plot')
xlabel('iteration number')
ylabel('deltaX value')
end

