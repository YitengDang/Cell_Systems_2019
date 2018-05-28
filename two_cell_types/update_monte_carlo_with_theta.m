function [p_out, theta_out, pe, cont] = update_monte_carlo_with_theta(N1, N2, p_in, theta_in, Con, Coff, K, f_mat, g_mat)
%%
% Updates cells with a Monte Carlo step
% array entries = cell types

% Required format:
% p = [p1; p2]
% theta = [theta_11 theta_12 theta_22];
% Con = [Con1 Con2]
% K = [K1 K2]
% f_mat = [f11 f12; f21 f22];
% g_mat = [g11 g12; g21 g22];

%p = [Non1(end)/N1; Non2(end)/N2];
%p_in = [0; 0];
%[px, py] = meshgrid(p_in, p_in);
%theta_in = f_mat.*(2*px-1).*(2*py-1);

%theta = [theta_11 theta_12 theta_22];

% check input formats
if ~(all(size(p_in)==[2 1]) && all(size(theta_in)==[2 2]) && all(size(Con)==[1 2])...
        && all(size(Coff)==[1 2]) && all(size(K)==[1 2]) && all(size(f_mat)==[2 2])...
        && all(size(g_mat)==[2 2]) )
    warning('Input has wrong size! Operation terminated');
    disp('Input has wrong size! Operation terminated');
    cont = 0;
    return
end

N = N1+N2;
N_12 = [N1; N2];
Non = round([N1; N2].*p_in);

% correct Con, Coff
Mcomm = f_mat>0;
Con = diag(Mcomm)'.*Con;
Coff = diag(Mcomm)'.*Coff;

%theta_mat = [theta(1) theta(2); theta(2) theta(3)]; 
xl = [N1/N; N2/N];

%---Calculate Pon->on, Poff->off---
%Xlm_on = diag(1./(2.*p))*(theta_mat + f_mat.*diag(2*p-1));
%Xlm_off = diag(1./(2.*(1-p)))*(-theta_mat + f_mat.*diag(2*p-1));
Xlm_on = diag(1./(2.*p_in))*(theta_in + repmat(2*p_in'-1, 2, 1));
Xlm_off = diag(1./(2.*(1-p_in)))*(-theta_in + repmat(2*p_in'-1, 2, 1));

S_on = repmat((Con-Coff)/2, 2, 1).*Xlm_on + repmat((Con+Coff)/2, 2, 1);
S_off = repmat((Con-Coff)/2, 2, 1).*Xlm_off + repmat((Con+Coff)/2, 2, 1);

muon = sum(diag(1./xl)*f_mat.*S_on, 2);
muoff = sum(diag(1./xl)*f_mat.*S_off, 2);

S2 = ((Con'-Coff').^2.*p_in.*(1-p_in));
sigmaon = sqrt(diag(1./xl)*g_mat*S2);
sigmaoff = sigmaon;

zon = (K' - Con' - muon)./sigmaon;
zoff = (K' - Coff' - muoff)./sigmaoff;
Poffoff = normcdf(zoff);
Ponon = 1-normcdf(zon);

%%
%---Update p1, p2---
yonon = binornd(Non, Ponon); % ON->ON
yoffoff = binornd(N_12-Non, Poffoff); % ON->ON

yoffon = (N_12 - Non) - yoffoff; % OFF -> ON
yonoff = Non - yonon; % ON->OFF

delta_p = (yoffon - yonoff)./N_12;

% update p
p_out = p_in + delta_p; % = [p1_out p2_out]
p_out = max(min(p_out, 1), 0);
%%
%---Calculate termination probability---
% exceptions if one of the probabilities is not defined (because there are no cells of
% a certain type)
if sum(isnan(Poffoff))>0 
    pe = Ponon.^N_12;
elseif sum(isnan(Ponon))>0
    pe = Poffoff.^N_12;
else
    pe = (Ponon.^p_in.*Poffoff.^(1-p_in)).^N_12;
end
pe = pe';

% continue?
if rand() < prod(pe)
    cont = 0; 
else
    cont = 1;
end
%% Calculate final theta 
% set NaN values to 1 for calculation of theta
Ponon(isnan(Ponon)) = 1;
Poffoff(isnan(Poffoff)) = 1;

theta_2 = theta_in; % temporary variable
theta_2(f_mat==0) = 0; % set theta values for which f^lm = 0 equal to zero.

%------v2: mean-field estimate for Delta 2---------------------------------
%delta_2 = 4*delta_p(1)*f_mat*delta_p(2);
%--------------------------------------------------------------------------
p_2 = repmat(p_in, 1, 2);
Ponon_2 = repmat(Ponon, 1, 2);
Poffoff_2 = repmat(Poffoff, 1, 2);
delta_p_2 = repmat(delta_p, 1, 2);

%---v3: exact solution for Delta_2-----------------------------------------
% Doesn't seem to work because possibly denom = 0
%{
try
    delta_2_nom = (-2.*f_mat.^2.*(-1 + 2.*p_2).*(Poffoff_2 - Ponon_2).*(-1 + p_2 + delta_p_2).*(p_2 + delta_p_2)...
        + f_mat.*delta_p_2.*(-1 + Ponon_2 + (Poffoff_2 - Ponon_2).*(p_2 + delta_p_2))...
        + 2.*f_mat.*(-2 + Poffoff_2 + Ponon_2).*(-1 + p_2 + delta_p_2).*(p_2 + delta_p_2).*theta_2...
        - delta_p_2.*(1 - Ponon_2 + p_2.*(-2 + Poffoff_2 + Ponon_2)...
        + (-2 + Poffoff_2 + Ponon_2).*delta_p_2).*theta_2);
    delta_2_denom = ((-1 + p_2).*p_2.*(2 + f_mat - Poffoff_2 - Ponon_2)...
        + (-1 + f_mat.*(-1 + 2.*p_2) + Poffoff_2 - p_2.*(-2 + Poffoff_2 + Ponon_2)).*delta_p_2...
        + f_mat.*delta_p_2.^2);
    delta_2 = delta_2_nom./delta_2_denom;
catch
    warning('delta_2');
end
%}
%----v3b: exact solution that takes into account p==0 and p==1 cases
try 
    delta_2_nom = (-2.*f_mat^2.*(-1 + 2.*p_2).*(Poffoff_2 - Ponon_2).*(-1 + p_2 + delta_p_2).*(p_2 + delta_p_2)...
        - delta_p_2.*(1 - Ponon_2 + p_2.*(-2 + Poffoff_2 + Ponon_2) + (-2 + Poffoff_2 + Ponon_2).*delta_p_2).*theta_in ...
        + f_mat.*(-p_2 + p_2^2 - 2.*delta_p_2 + 2.*p_2.*delta_p_2 + p_2.*Poffoff_2.*delta_p_2 + Ponon_2.*delta_p_2 - p_2.*Ponon_2.*delta_p_2 + delta_p_2^2 ...
        + Poffoff_2.*delta_p_2^2 - Ponon_2.*delta_p_2^2 + 2.*(-2 + Poffoff_2 + Ponon_2).*(-1 + p_2 + delta_p_2).*(p_2 + delta_p_2).*theta_in)...
        - (-1 + p_2 + delta_p_2).*(f_mat.*(p_2 + 2.*delta_p_2 - Ponon_2.*delta_p_2) - ...
        (-1 + Ponon_2).*delta_p_2.*theta_in).*(p_2==0) +...
        (-1 + Poffoff_2).*(p_2 + delta_p_2).*(-f_mat.*delta_p_2 + 2.*f_mat^2.*(-1 + 2.*p_2).*(-1 + p_2 + delta_p_2)...
        + delta_p_2.*theta_in - 2.*f_mat.*(-1 + p_2 + delta_p_2).*theta_in).*(p_2==1));
    
    delta_2_denom = (-1 + p_2).*p_2.*(2 + f_mat - Poffoff_2 - Ponon_2)...
        + (-1 + f_mat.*(-1 + 2.*p_2) + Poffoff_2 - p_2.*(-2 + Poffoff_2 + Ponon_2)).*delta_p_2 + f_mat.*delta_p_2^2 ...
        +  p_2.*(-1 + Ponon_2).*(-1 + p_2 + delta_p_2).*(p_2==0)...
        + (-1 + p_2).*(-1 + Poffoff_2).*(p_2 + delta_p_2).*(p_2==1);
    
    delta_2 = delta_2_nom./delta_2_denom;
    
    delta_2_alt = 4*delta_p(1)*f_mat*delta_p(2);
    delta_2(isnan(delta_2)) = delta_2_alt(isnan(delta_2)); % temporary solution: set Delta 2 = mean-field for undefined
catch
    warning('Delta 2 calculation wrong');
end


%--------------------------------------------------------------------------

dtheta =...
    -(theta_2.*f_mat + diag(2*p_in-1)*f_mat).*(1-repmat(Ponon, 1, 2))... 
    -(theta_2.*f_mat - diag(2*p_in-1)*f_mat).*(1-repmat(Poffoff, 1, 2))...
    -(theta_2.*f_mat + f_mat*diag(2*p_in-1)).*(1-repmat(Ponon', 2, 1))...
    -(theta_2.*f_mat - f_mat*diag(2*p_in-1)).*(1-repmat(Poffoff', 2, 1))...
    + delta_2;

theta_out = theta_in;
idx = f_mat>0;
theta_out(idx) = theta_in(idx) + dtheta(idx)./f_mat(idx); %theta is normalized by f_ij

% correct range
theta_out = max(min(theta_out, ones(2,2)), -ones(2,2)); 
%% Previous code

%dtheta_11 = -1/xl(1)*(theta_mat(1,1) + (2*p(1)-1)*f_mat(1,1))*(1-Ponon(1)) ... 
%    -1/xl(1)*(theta_mat(1,1) - (2*p(1)-1)*f_mat(1,1))*(1-Poffoff(1))...
%    -1/xl(1)*(theta_mat(1,1) + (2*p(1)-1)*f_mat(1,1))*(1-Ponon(1))...
%    -1/xl(1)*(theta_mat(1,1) + (2*p(1)-1)*f_mat(1,1))*(1-Poffoff(1))...
%    + 4*f_mat(1,1)*delta_p(1)*delta_p(2);

%dtheta_12 = -1/xl(1)*(theta_mat(1,2) + (2*p(1)-1)*f_mat(1,2))*(1-Ponon(1)) ... 
%    -1/xl(1)*(theta_mat(1,2) - (2*p(1)-1)*f_mat(1,2))*(1-Poffoff(1))...
%    -1/xl(2)*(theta_mat(1,2) + (2*p(2)-1)*f_mat(1,2))*(1-Ponon(2))...
%    -1/xl(2)*(theta_mat(1,2) + (2*p(2)-1)*f_mat(1,2))*(1-Poffoff(2))...
%    + 4*f_mat(1,2)*delta_p(1)*delta_p(2);

%dtheta =...
%    -diag(1./xl)*(theta_in.*f_mat + diag(2*p_in-1)*f_mat).*(1-Ponon_mat_1)... 
%    -diag(1./xl)*(theta_in.*f_mat - diag(2*p_in-1)*f_mat).*(1-Poffoff_mat_1)...
%    -((theta_in.*f_mat + f_mat*diag(2*p_in-1)).*(1-Ponon_mat_2))*diag(1./xl)...
%    -((theta_in.*f_mat - f_mat*diag(2*p_in-1)).*(1-Poffoff_mat_2))*diag(1./xl)...
%    + 4*delta_p(1)*f_mat*delta_p(2);
%dtheta =...
%    -diag(1./xl)*(theta_2.*f_mat + diag(2*p_in-1)*f_mat).*(1-Ponon_mat_1)... 
%    -diag(1./xl)*(theta_2.*f_mat - diag(2*p_in-1)*f_mat).*(1-Poffoff_mat_1)...
%    -((theta_2.*f_mat + f_mat*diag(2*p_in-1)).*(1-Ponon_mat_2))*diag(1./xl)...
%    -((theta_2.*f_mat - f_mat*diag(2*p_in-1)).*(1-Poffoff_mat_2))*diag(1./xl)...
%    + 4*delta_p(1)*f_mat*delta_p(2);

% L=1
% muon = fN*(Con*p + 1 - p + (Con-1)*(1-p)*I)
% muoff = fN*(Con*p + 1 - p - (Con-1)*p*I);
