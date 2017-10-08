%% Calculate the derivative
% This function gives will use either forward difference or central
% differencing depending on user input to evaluate the first derivative of
% objective function, contraint function and augmented lagrangian function

function [L_deriv, f_deriv, c_deriv] =  FirstDerivativeCalculations(f,c,L,epi,x_parameter,r0,v0,lambda,rho,scaleFactor,FDS,CDS)

% Copy the orginal for derivative calculations later on
f0 = f;
c0 = c;
L0 = L;
x0 = x_parameter;

% Allocate matrix 
f_deriv = zeros(size(x_parameter))';
c_deriv = zeros(length(lambda), size(x_parameter,2));
L_deriv = zeros(size(x_parameter))';

%% Forward differencing scheme
% If user chose forward differencing scheme
if FDS == 1
    for i = 1:length(x_parameter)
        
        % Extract orginal states
        x_parameter = x0;
        % Add forward perturbation to the particular state that "i" is on
        x_parameter(i) = x_parameter(i) + epi;
        
        % Evaluate the new objective function, constraint
        % function and augmented lagarangian function with
        % perturbed state
        [LPerturbed, cPerturbed, fPerturbed] = BurnDynamicModel(x_parameter,r0,v0,scaleFactor,lambda,rho);
        
        % Apply Forward Differncing Scheme
        f_deriv(i) = (fPerturbed - f0)/epi;
        c_deriv(:,i) = (cPerturbed - c0)./epi;
        L_deriv(i) = (LPerturbed - L0)/epi;
    end
end
%% Central differncing scheme
% If the user chose central differencing scheme
if CDS == 1
    for i = 1:length(x_parameter)
        
        % Extract orginal states and put into forward and backward
        x_Forward = x0;
        x_Backward = x0;
        % Add forward perturbation to the particular state that "i" is on
        x_Forward(i) = x_parameter(i) + epi;
        % Add backward perturbation to the particular state that "i" is on
        x_Backward(i) = x_parameter(i) - epi;
        
        % Evaluate the new objective function, constraint
        % function and augmented lagarangian function with each
        % perturbed state. Therefore, a forward and a backward one are found
        [LPerturbed_F, cPerturbed_F, fPerturbed_F] = BurnDynamicModel(x_Forward,r0,v0,scaleFactor,lambda,rho);
        [LPerturbed_B, cPerturbed_B, fPerturbed_B] = BurnDynamicModel(x_Backward,r0,v0,scaleFactor,lambda,rho);
        
        % Apply Forward Differncing Scheme
        f_deriv(i) = (fPerturbed_F - fPerturbed_B)/(2*epi);
        c_deriv(:,i) = (cPerturbed_F - cPerturbed_B)./(2*epi);
        L_deriv(i) = (LPerturbed_F - LPerturbed_B)/(2*epi);
        
    end
end
end