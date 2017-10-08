%% Hessian Approximation using second order differencing
% This function implments aprroximation of Hessian by forward differencing
% the gradient. General formula obtained from lecture slide (Numerical
% Derivatives)
%
% Author: Kuan Chun Hwang

function H = SecondOrderHessianApprox(L0,epi,x_parameter,r0,v0,lambda,rho,scaleFactor)

% Store varaibles into i, j and ij where ij takes account into the combined
% effect of propergation
x0 = x_parameter;
xei = x_parameter;
xej = x_parameter;
xei_ej = x_parameter;

% Initialize the Hessian matrix which has the size of 8x8
H = zeros(8,8);

for i = 1:length(x_parameter)
    
    for j = 1:length(x_parameter)
        
        xei(i) = x_parameter(i) + epi;
        [L_ei,c_ei,f_ei] = BurnDynamicModel(xei,r0,v0,scaleFactor,lambda,rho);
        
        xei_ej(j) = xei(j) + epi;
        [L_eij,c_eij,f_eij] = BurnDynamicModel(xei_ej,r0,v0,scaleFactor,lambda,rho);
        
        xej(j) = x_parameter(j) + epi;
        [L_ej,c_ej,f_ej] = BurnDynamicModel(xej,r0,v0,scaleFactor,lambda,rho);
        
        % Compute new Hessian matrix
        H(j,i) = (L_eij - L_ei - L_ej + L0)/(epi^2);
        
        % Restore to original
        xej = x0;
        xei_ej(j) = xei_ej(j) - epi;
    end
    
    % Restore to original
    xei = x0;
end