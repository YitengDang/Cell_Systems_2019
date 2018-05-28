function [p_out, terminate] = update_monte_carlo(N1, N2, p, Con, Coff, K, f_mat, g_mat)
%%
% Updates cells with a Monte Carlo step
% array entries = cell types
% muon = fN*(Con*p + 1 - p + (Con-1)*(1-p)*I)
% muoff = fN*(Con*p + 1 - p - (Con-1)*p*I);
%p = p_in; %[p1 p2];
N = N1+N2;
N_12 = [N1 N2];
Non = round([N1 N2].*p);

%---Calculate Pon->on, Poff->off---
%muon = [4 5.2]; %[muon1 muon2]; 
%muoff = [2.3 2.1]; %[muoff1 muoff2];
%sigmaon = [1.2 1.3]; %[sigmaon1 sigmaon2];
%sigmaoff = [0.5 0.6]; %[sigmaoff1 sigmaoff2];

S = ((Con-Coff)/2.*(2*p-1) + (Con+Coff)/2)';
muon = ([N1/N; N2/N].*f_mat*S)';
muoff = muon;

S2 = ((Con-Coff).^2.*p.*(1-p))';
sigmaon = sqrt([N1/N; N2/N].*g_mat*S2)';
sigmaoff = sigmaon;

zon = (K - Con - muon)./sigmaon;
zoff = (K - 1 - muoff)./sigmaoff;
Poffoff = normcdf(zoff);
Ponon = 1-normcdf(zon);
%%
%---Update p1, p2---
yonon = binornd(Non, Ponon); % ON->ON
yoffoff = binornd(N_12-Non, Poffoff); % ON->ON

yoffon = (N_12-Non) - yoffoff; % OFF -> ON
yonoff = Non - yonon; % ON->OFF

delta_p = (yoffon - yonoff)./N_12;

p_out = p + delta_p; % = [p1_out p2_out]
p_out = max(min(p_out, 1), 0);

%---Calculate termination probability---
pe = (Ponon.^p.*Poffoff.^(1-p)).^N_12;

if all(rand(2,1) < pe)
    terminate = 1;
else
    terminate = 0;
end

