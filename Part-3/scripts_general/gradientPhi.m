%% Calculate the gradient of phi
%
% Author Kuan Chun Hwang

function gradPhi = gradientPhi( phihdl,alpha,phi,pert)

    for i=1:length(alpha)
        
        % Store orginal aplha
        alphaPert=alpha;
        
        % Add forward perturnation
        alphaPert(i) = alphaPert(i) + pert;     
        
        % Find gradient
        gradPhi = (phihdl(alphaPert) - phi)/pert;
        
    end
end