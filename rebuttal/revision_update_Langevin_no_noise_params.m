function [p_new, I_new, Peq, p_lim] = revision_update_Langevin_no_noise_params(pt, It, N, Con, K, fN, gN, alpha, M_int)
    % Estimate delta from <Delta p>
    if nargin<9
        M_int = 1;
    end
    
    delh_delp = @(p,I) M_int*(Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;
    delh_delI = @(p,I) M_int*2*fN*(Con - 1)*p.*(p - 1);

    muon = fN*(Con*pt + 1 - pt + (Con-1)*(1-pt)*It);
    muoff = fN*(Con*pt + 1 - pt - (Con-1)*pt*It);
    kappa = (Con-1)*sqrt(gN*pt*(1-pt));
    sigmaon = sqrt(kappa^2+alpha^2);
    sigmaoff = sqrt(kappa^2+alpha^2);
    
    zon = (K - Con - muon)/sigmaon;
    zoff = (K - 1 - muoff)/sigmaoff;
    if M_int==1
        Poffoff = normcdf(zoff);
        Ponon = 1-normcdf(zon);

        Peq = (Ponon.^pt .* Poffoff.^(1-pt)).^N;
    
    elseif M_int==-1
        Poffoff = 1-normcdf(zoff);
        Ponon = normcdf(zon);

        Peq = (Ponon.^pt .* Poffoff.^(1-pt)).^N;
    end
    
    delp = (1-pt)*(1-Poffoff) - pt*(1-Ponon);
    
    % match noise with std(Delta p)
    delta = -delp/delh_delp(pt, It);
    noise = sqrt(pt*(1-Poffoff)*Poffoff + pt*(1-Ponon)*Ponon)/sqrt(N);
    
    pnoise = normrnd(0, noise);
    Inoise = pnoise*delh_delI(pt, It)/delh_delp(pt, It);
    p_new = pt + delp + pnoise;
    %I_new = It + delta*delh_delI(pt, It)/delh_delp(pt, It)) + Inoise;
    %p_new = pt - delh_delp(pt, It)*delta + pnoise;
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