function [p_new, I_new, Peq, p_lim] = EOM_update_Langevin(pt, It, N, Con, K, fN, gN, delta, sigma, alpha)
    % New version: modified drift for I     
    delh_delp = @(p,I) (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;
    delh_delI = @(p,I) 2*fN*(Con - 1)*p.*(p - 1);

    muon = fN*(Con*pt + 1 - pt + (Con-1)*(1-pt)*It);
    muoff = fN*(Con*pt + 1 - pt - (Con-1)*pt*It);
    kappa = sqrt((Con-1)^2*(gN)*pt*(1-pt));
    sigmaon = sqrt(kappa^2+alpha^2);
    sigmaoff = sqrt(kappa^2+alpha^2);

    zon = (K - Con - muon)/sigmaon;
    zoff = (K - 1 - muoff)/sigmaoff;
    Poffoff = normcdf(zoff);
    Ponon = 1-normcdf(zon);

    Peq = (Ponon^pt*Poffoff^(1-pt))^N;
    
    pnoise = normrnd(0,sigma);
    Inoise = normrnd(0,sigma)*delh_delI(pt, It)/delh_delp(pt, It);
    p_new = pt - delh_delp(pt, It)*delta + pnoise;
    I_new = It - delh_delI(pt, It)*delta + Inoise;
    p_lim = 0;
    
    if p_new <= 0
       p_new = 0;
       I_new = 0;
       p_lim = -1;
       %Peq = 1;
    elseif p_new >= 1
       p_new = 1;
       I_new = 0;
       p_lim = 1;
       %Peq = 1;
    end
end