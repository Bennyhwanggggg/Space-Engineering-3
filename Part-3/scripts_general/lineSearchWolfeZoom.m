%% function 'lineSearchWolfeZoom'
%

function alphaOpt = lineSearchWolfeZoom(alpha_lo, alpha_hi,...
        phihdl, phidhdl, c1, c2)
    

    % Maximum number of iterations
    jMax = 100;

    alpha2 = zeros(jMax, 1);
    phi2 = zeros(jMax, 1);
    phid2 = zeros(jMax, 1);

    j = 1;

    alphaOpt = [];
    phi0 = phihdl(0);
    phid0 = phidhdl(0);
    % selection phase
    while j < jMax
        alpha2(j) = lineSearchWolfeInterpolate(alpha_lo, alpha_hi, ...
            phihdl(alpha_lo), phihdl(alpha_hi), ...
            phidhdl(alpha_lo), phidhdl(alpha_hi));
        phi2(j) = phihdl(alpha2(j)); % evaluate phi(j)   

        if ( (phi2(j)>phi0+c1*alpha2(j)*phid0) ...
                || (phi2(j)>=phihdl(alpha_lo)) )
            % sufficient decrease condition violated
            alpha_hi = alpha2(j);
        else
            phid2(j) = phidhdl(alpha2(j)); % evaluate phi(j)   
            if (abs(phid2(j)) <= -c2*phid0)
                % solution found
                alphaOpt = alpha2(j);
                break;
            end

            if (phid2(j)*(alpha_hi - alpha_lo) >= 0)
                alpha_hi = alpha_lo;
            end
            alpha_lo = alpha2(j);
        end
        j = j + 1;
    end

    alpha2 = alpha2(1:j, 1); % remove entries
    
    if isempty(alphaOpt)
        
        alphaOpt = alpha2(end);
    end
 
end
