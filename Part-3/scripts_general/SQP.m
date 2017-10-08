%% Sequential Qudartic Programing (SQP)
%
% Author Kuan Chun Hwang

function [x_parameter,L_data,c_data,f_data,alpha_data,error_data,x_data,count_data] = SQP(x_parameter,r0,v0,scaleFactor,FDS, CDS, SOD, BFGS, amer, bmer)
% Intitialze SQP variables
H       = eye(8,8);     % Initial Hessian guess matrix as an 8x8 identity matrix
tol     = 0.0001;       % Set tolerance value
epi    = 0.000001;      % Set perturbation to a small value
lambda  = zeros(5,1);   % Initial Lagrangian Multipliers
rho     = 10000;        % Set augmented Lagrangian penalty term
alpha_initial = 1;      % Initialize back tracking step
error   = 1000;         % Initialize error
count   = 1;            % Initialize count

while error > tol
    
    % Calculate the objective function, contraint function and augmented
    % lagrangian
    [L,c,f] = BurnDynamicModel(x_parameter,r0,v0,scaleFactor,lambda,rho);
    
    % Store the objective, constraint and Lagrangian
    L_data(:,count) = L;
    f_data(:,count) = f;
    c_data(:,count) = c;
    
    % Compute first derivative
    [L_deriv, f_deriv, c_deriv] =  FirstDerivativeCalculations(f,c,L,epi,x_parameter,r0,v0,lambda,rho,scaleFactor,FDS,CDS);
    
    % Check if Hessian is positive definite
    if(min(eig(H)) < 0)
        
        % If not, set it back to an identitiy matrix of 8x8
        H = eye(8,8);
        
    end
    
    % Find KKT system LHS
    G = c_deriv;
    KKT_LHS = [H, G'; G, zeros(5,5)'];
    % Find KKT system RHS
    g = f_deriv';
    KKT_RHS = [g; c];
    % Solve the KKT system
    negP_lambdaplus = KKT_LHS\KKT_RHS;
    % Extract search direction
    P = -negP_lambdaplus(1:length(x_parameter));
    
    % Use the line serach algorithm to find step size alpha and store step
    % size
    alpha = lineSearchWolfeAlgorithm(x_parameter, r0, v0, P, scaleFactor, lambda, rho, epi, alpha_initial);
    alpha_data(:,count) = alpha;
    
    % Find new x
    x_parameterNew = x_parameter + alpha*P;
    
    % Use either merit function or not depending on user input
    if amer == 1
        % Calculate new lambda value for merit function
        lambda = ((G*G')\G)*g;
    elseif bmer == 1
        % lambda for second order differencing
        lambda = negP_lambdaplus(9:end);
    end
    
    % Chose the Hessian approximation method depending on user input
    if SOD == 1
        % Use second order differencing to approximate Hessian matrix
        H = SecondOrderHessianApprox(L,epi,x_parameter,r0,v0,lambda,rho,scaleFactor);
    elseif BFGS == 1 && amer == 1 && CDS ==1
        % Use Quasi Newtom -> BFGS method to approximate Hessian matrix
        L_derivNew = g-G'*lambda;
        H = QuasiNewtonBFGS(H,x_parameterNew-x_parameter,L_derivNew-L_deriv');
    end
    
    % Compute error and store error in an array
    error = abs(max(x_parameterNew-x_parameter));
    error_data(:,count) = error;
    
    % Store x data to help keep track of optimisation progress and update the parameters
    x_data(:,count) = x_parameter;
    x_parameter = x_parameterNew;
    
    % Update count
    count_data(:,count) = count;
    count = count +1;
end
end