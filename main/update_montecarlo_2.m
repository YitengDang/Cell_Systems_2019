function [theta_new, p_new, pe] = update_montecarlo_2(theta, p, N, Con, K, fN, gN, alpha)
% v2: test an alternative formula for p
    phi = @(x) exp(-0.5*x^2)/sqrt(2*pi);
    if p <= 0 || p >= 1
        I = 0;
    else
        I = (theta - (2*p-1)^2*fN)/4/p/(1-p)/fN;
    end
    muon = fN*(Con*p + 1 - p + (Con-1)*(1-p)*I);
    muoff = fN*(Con*p + 1 - p - (Con-1)*p*I);
    kappa = sqrt((Con-1)^2*(gN)*p*(1-p));
    sigmaon = sqrt(kappa^2+alpha^2);
    sigmaoff = sqrt(kappa^2+alpha^2);

    zon = (K - Con - muon)/sigmaon;
    zoff = (K - 1 - muoff)/sigmaoff;
    Poffoff = normcdf(zoff);
    Ponon = 1-normcdf(zon);
    
    pe = (Ponon^p*Poffoff^(1-p))^N;

    % Calculate the average fraction of cells that change state and the
    % variance
%     p_minus_mean = p*(1-Ponon);
%     p_plus_mean = (1-p)*(1-Poffoff);
%     
%     var_p_plus = p_plus_mean*Poffoff/N;
%     var_p_minus = p_minus_mean*Ponon/N;
%     
%     if var_p_plus > 0
%         p_plus = normrnd(p_plus_mean, sqrt(var_p_plus));
%     else
%         p_plus = p_plus_mean;
%     end
%     
%     if var_p_minus > 0
%         p_minus = normrnd(p_minus_mean, sqrt(var_p_minus));
%     else
%         p_minus = p_minus_mean;
%     end
    
    % Calculate number of cells that turn ON and OFF respectively
    p_minus = binornd(round(N*p), 1-Ponon)/N;
    p_plus = binornd(round(N*(1-p)), 1-Poffoff)/N;
    
    % Calculate the average and variance of dtheta. Take care that the
    % function has a problem when P(off->off) and/or P(on->on) are close to
    % unity.
    naux = round(p_minus*N);
    if isreal(sigmaon) && sigmaon > 0 && naux > 0
        zon = (K - Con - muon)/kappa*ones(naux,1);
        Y_minus = mean(muon + ...
            kappa*trandn((fN-muon)/kappa*ones(naux,1),zon));
    else
        Y_minus = 0;
    end
    
    naux = round(p_plus*N);
    if isreal(sigmaoff) && sigmaoff > 0 && naux > 0
        zoff = (K - 1 - muoff)/kappa*ones(naux,1);
        Y_plus = mean(muoff + ...
            kappa*trandn(zoff, (Con*fN-muoff)/kappa*ones(naux,1)));
    else
        Y_plus = 0;
    end
    
    % Evolve the state one time step
    dp = -(1-Ponon)*p+(1-Poffoff)*(1-p);
    %dp = p_plus - p_minus;
    dtheta = 8/(Con-1)*(p_plus*Y_plus - p_minus*Y_minus) - ...
        4*(Con+1)/(Con-1)*fN*dp + 4*fN*dp^2;

    theta_new = theta + dtheta;
    
    theta_new = min(theta_new, fN);
    if theta_new == fN
        p_new = round(p);
    else
        p_new = p + dp;
        p_new = round(N*p_new)/N; % To make it a possible value
        p_new = max(0,min(p_new,1));
    end