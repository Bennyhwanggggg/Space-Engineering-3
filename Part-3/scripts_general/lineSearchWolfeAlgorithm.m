%% Line Search Algorithm for the Wolfe conditions
% Implements the line search algorithm using a zoom function such that both
% Wolfe conditions are satisfied
%
% Auhtor: Kuan Chun Hwang

function alphaOpt = lineSearchWolfeAlgorithm(x_parameter, r0, v0, P, scaleFactor, lambda, rho, pert, alpha_initial)

% Initialize the Wolfe conditions variables and handle functions
% Strong Wolfe condition satisfies 0<c1<c2<1
c1 = 0.0001;
c2 = 0.8;
% Initalize alpha such that alpha0 <- 0  and alphaMax > 0 and alpha1 =>
% 0~alphaMax. In this case, the index number 1 is used as a 0.
alpha(1) = 0;
alphaMax = 1;

% Initialize function handle which finds the lagrangian function and its derivative
phihdl = @(alpha) BurnDynamicModel(x_parameter+alpha*P,r0,v0,scaleFactor,lambda,rho);
phidhdl = @(alpha) gradientPhi( phihdl,alpha,phihdl(alpha),pert);

% Preallocate phi and dhpi
phi = zeros(100,1);
dphi = zeros(100,1);
% To peform iteration, we are required to find phi(1) and its gradient
phi(1) = phihdl(alpha(1));
dphi(1) = phidhdl(alpha(1));
% Since 1 is taken as 0, the second index is taken as the initial aplha 1
% value which is within the range of 0~alphaMax
alpha(2) = alpha_initial;

% Start iteration at the second index since we are required to search for
% phi(alpha(i-1)) according to pseudo code but we do not have anything at alpha(0)
i = 2;

% Set maximum iteration number
iter_max = 100;

% Iteration process
while i < iter_max
    
    % Compute phi at alpha(i)
    phi(i) = phihdl(alpha(i));
    if phi(i) > phi(1) + c1*alpha(i)*dphi(1) || (phi(i) >= phi(i-1) && i>2)
        
        % Use zoom function provided by tutors to find alpha and stop
        alphaOpt = lineSearchWolfeZoom(alpha(i-1),alpha(i),...
            phihdl, phidhdl, c1, c2);
        % Stop when alpha is found
        return;
    end
    
    % Compute dphi at alpha(i)
    dphi(i) = phidhdl(alpha(i));
    if abs(dphi(i)) <= -c2*dphi(1)
        % Set alphaOPt as the alpha(i) and stop
        alphaOpt = alpha(i);
        % Stop when alpha is found
        return;
    end
    if dphi(i) >= 0
        % Use zoom function provided by tutors to find alpha and stop
        alphaOpt = lineSearchWolfeZoom(alpha(i-1),alpha(i),...
            phihdl,phidhdl,c1,c2);
        % Stop when alpha is found
        return;
    end
    % Chose the next alpha which has to be within the range of alpha(i) and alphamax.
    alpha(i+1) = (alpha(i)+alphaMax)/2;
    i=i+1;
end

% After 100 iteration, if no alphaOpt is found, use the alpha from last
% iteration
alphaOpt = alpha(i);
end
