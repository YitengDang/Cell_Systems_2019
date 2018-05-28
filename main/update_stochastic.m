function [theta_new, p_new, pe] = update_stochastic(theta, p, N, Con, K, fN, gN, alpha)
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
    ponoff = p*(1-Ponon);
    poffon = (1-p)*(1-Poffoff);
    dp = poffon-ponoff;
    vardp = (poffon*Poffoff + ponoff*Ponon)/N;
    
    if vardp < 0
        vardp = 0;
    end
    
    dtheta = 8/(Con-1)*(poffon*muoff - ponoff*muon + ...
        sigmaoff*(1-p)*phi(zoff)+ sigmaon*p*phi(zon)) ...
        - 4*fN*(Con+1)/(Con-1)*dp + 4*fN*(dp^2 + vardp);

    vartheta = 64/(Con-1)^2*((1-p)^2*sigmaoff^2* ...
        ((1-Poffoff)^2+zoff*phi(zoff)*(1-Poffoff) - (phi(zoff))^2) ...
        + p^2*sigmaon^2* ...
        ((1-Ponon)^2-zon*phi(zon)*(1-Ponon) - (phi(zon))^2)) + ...
        16*fN^2*(4*dp^2*vardp + 2*vardp^2 + (Con+1)^2/(Con-1)^2*vardp);

    if vartheta < 0
        vartheta = 0;
    end
    
    pe = (Ponon^p*Poffoff^(1-p))^N;

    p_new = round(N*(p + normrnd(dp, sqrt(vardp))))/N;
    p_new = max(0,min(p_new,1));
    theta_new = theta + normrnd(dtheta, sqrt(vartheta));
    theta_new = min(theta_new, fN);
    if theta_new == fN
        p_new = round(p_new);
    end
    if p_new == 0 || p_new == 1
        theta_new = fN;
    end