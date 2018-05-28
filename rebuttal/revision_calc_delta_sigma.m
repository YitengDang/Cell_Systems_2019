function [delta, sigma] = revision_calc_delta_sigma(N, K, Con, fN, gN)
% Estimate delta and sigma in Langevin equation
    % Calculate as for loop
    delta = 0;
    sigma = 0;
    I = 0;

    for n=0:N
        p = n/N;
        % Gradient term
        delh_delp = (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;

        % Calculate Ponon, Poffoff => delp
        muon = fN*(Con*p + 1 - p + (Con-1)*(1-p)*I);
        muoff = fN*(Con*p + 1 - p - (Con-1)*p*I);
        kappa = (Con-1)*sqrt(gN*p*(1-p));
        sigmaon = kappa; %sqrt(kappa^2+alpha^2);
        sigmaoff = kappa; %sqrt(kappa^2+alpha^2);

        zon = (K - Con - muon)/sigmaon;
        zoff = (K - 1 - muoff)/sigmaoff;
        Poffoff = normcdf(zoff);
        Ponon = 1-normcdf(zon);

        delp = (1-p)*(1-Poffoff) - p*(1-Ponon);
        prob = nchoosek(N, n)/2^N;

        % Update delta, sigma
        delta = delta - prob*delp/delh_delp;
        sigma = sigma + prob*sqrt(p*(1-Poffoff)*Poffoff + p*(1-Ponon)*Ponon)/sqrt(N);
    end
end