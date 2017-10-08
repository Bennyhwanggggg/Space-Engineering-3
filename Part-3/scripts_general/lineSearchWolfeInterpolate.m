%% function line search wolfe interpolation
%
% 

function alphaInterp = lineSearchWolfeInterpolate(alpha0, alpha1,...
    phi0, phi1, phid0, phid1)

    % 3rd order polynomial interpolation
    d1 = phid0 + phid1 - 3*(phi0 - phi1)/(alpha0 - alpha1);
    d2 = sign(alpha1 - alpha0)*sqrt(d1^2 - phid0*phid1);
    alphaInterp = alpha1 - (alpha1 - alpha0)*(phid1 + d2 - d1)...
        /(phid1 - phid0 + 2*d2); % minimum of the 3rd order model
end