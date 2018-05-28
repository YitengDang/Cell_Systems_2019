function [theta_new, p_new, pe] = update_stochastic_cov(theta, p, N, a0, Rcell,...
    Con, K, fN, gN, alpha)
    phi = @(x) exp(-0.5*x^2)/sqrt(2*pi);
    if p <= 0 || p >= 1
        I = 0;
    else
        I = (theta - (2*p-1)^2*fN)/4/p/(1-p)/fN;
    end
    muon = fN*(Con*p + 1 - p + (Con-1)*(1-p)*I);
    muoff = fN*(Con*p + 1 - p - (Con-1)*p*I);
    kappa = (Con-1)^2*(gN)*p*(1-p);
    sigmaon = sqrt(kappa+alpha^2);
    sigmaoff = sqrt(kappa+alpha^2);
    if ~isreal(sigmaon) || sigmaon == 0
        sigmaon = 0;
        if K - Con - muon > 0
            zon = inf;
        else
            zon = -inf;
        end
    else
        zon = (K - Con - muon)/sigmaon;
    end
    
    if ~isreal(sigmaoff) || sigmaoff == 0
        sigmaoff = 0;
        if K - 1 - muoff > 0
            zoff = inf;
        else
            zoff = -inf;
        end
    else
        zoff = (K - 1 - muoff)/sigmaoff;
    end
        
    Poffoff = normcdf(zoff);
    Ponon = 1-normcdf(zon);
    p_minus = p*(1-Ponon);
    p_plus = (1-p)*(1-Poffoff);
    
    dp = p_plus - p_minus;
    
    var_p_plus = p_plus*Poffoff/N;
    var_p_minus = p_minus*Ponon/N;
    vardp = var_p_plus + var_p_minus;
    
    if vardp < 0
        vardp = 0;
    end
    
    if 1 - Poffoff < 1e-16
        y_plus = muoff;
        var_y_plus = sigmaoff^2;
    else
        y_plus = muoff + sigmaoff*phi(zoff)/(1 - Poffoff);
        if sigmaoff < 1e-16
            var_y_plus = 0;
        else
            var_y_plus = sigmaoff^2*(1 + zoff*phi(zoff)/(1 - Poffoff) ...
                - (phi(zoff)/(1 - Poffoff))^2);
        end
    end
    
    if 1 - Ponon < 1e-16
        y_minus = muon;
        var_y_minus = sigmaon^2;
    else
        y_minus = muon - sigmaon*phi(zon)/(1 - Ponon);
        if sigmaon < 1e-16
            var_y_minus = 0;
        else
            var_y_minus = sigmaon^2*(1 - zon*phi(zon)/(1 - Ponon) ...
                - (phi(zon)/(1 - Ponon))^2);
        end
    end
    
    sigma_minus = p_minus^2*var_y_minus + ...
        y_minus^2*var_p_minus + var_y_minus*var_p_minus;
    
    sigma_plus = p_plus^2*var_y_plus + ...
        y_plus^2*var_p_plus + var_y_plus*var_p_plus;
    
    dtheta = 8/(Con-1)*(p_plus*y_plus - p_minus*y_minus) ...
        + 4*fN*(dp^2+vardp - dp*(Con+1)/(Con-1));
    
    vardtheta = 64/(Con-1)^2*(sigma_plus + sigma_minus) + ...
        16*fN^2*vardp*((Con+1)^2/(Con-1)^2 + 4*dp^2 + 2*vardp);
    
    covdpdtheta = 8/(Con-1)*(y_plus*(var_p_plus + p_plus*dp) ...
        + y_minus*(var_p_minus - p_minus*dp)) ...
        + 4*fN*(dp^3 + 3*dp*vardp - ...
        (Con+1)/(Con-1)*(vardp + dp^2)) - dp*dtheta;
    

    if vardtheta < 0
        vardtheta = 0;
    end
    
    cov_d = [vardp covdpdtheta; covdpdtheta vardtheta];
    mean_d = [dp; dtheta];
    
    pe = (Ponon^p*Poffoff^(1-p))^N;

    [c_aux, r] = chol(cov_d, 'lower');
    
    if r>0
        dx = [dp dtheta];
    else
        dx = mean_d + c_aux*polar_marsaglia;
    end
    p_new = round(N*(p + dx(1)))/N;
    theta_new = theta + dx(2);
    p_new = max(0,min(p_new,1));
    
    % Calculate maximum theta given a p
    theta_max = max_theta(p_new, a0, Rcell, N, fN);
    theta_new = min(theta_new, theta_max);
                