function [pe, p_lim] = revision_update_Peq(p, I, N, Con, K, fN, gN, alpha)
% calculate the equilibrium probability Pe
    pe = 1; % default
    if p <= 0
        %I = 0;
        p_lim = -1; % keeps track of whether the system has reached p=0 or p=1
    elseif p >= 1
        p_lim = 1;
    else
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
        p_lim = 0;
    end
end