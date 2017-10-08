%% Quasi Newton BFGS method to approximate Hessian
%
% Author: Kuan Chun Hwang

function Hnew = QuasiNewtonBFGS(H,sk,yk)

        % From formula, sk = delta X and yk = delta L
        Hnew = H - (H*(sk*sk')*H)/(sk'*H*sk) + (yk*yk')/(yk'*sk);
        
end