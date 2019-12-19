clear all
close all
clc
%% Setup
% parameters
network = 15;
gz = 15;
N = gz^2;
a0 = 0.5;
rcell = 0.2;
lambda = [1 1.2];

% get pos, dist
mcsteps = 0;
[~, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);
[fN, fnn, ~] = calc_fN(a0, rcell, gz, dist, lambda);

% mean-field 
pw = [1/8 1/8];
Rcell = rcell*a0;
f_MF = (4*pi/a0^2/sqrt(3)).*exp(Rcell./lambda).*sinh(Rcell./lambda).*exp(-3*a0./(2*lambda));

L = 1000;
N_inf = 0;

% save figure folder
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\robustness_analytical';
%%
showlines = 0;
h = plot_boundaries_K_Con(fN, fnn, pw, L, showlines);
[Q_value] = calc_Q_value(L, fN, fnn, pw, N_inf, f_MF)

% Calculate projected areas
A_frac = zeros(2);
A_frac(1,1) = calc_area(L, fN, fnn, pw, 'self', 1, 0, f_MF);
A_frac(1,2) = calc_area(L, fN, fnn, pw, 'mutual', 2, 0, f_MF);
A_frac(2,1) = calc_area(L, fN, fnn, pw, 'mutual', 1, 0, f_MF);
disp(A_frac)
%{
% Geometric check
C11_max = 459.9367;
C11_min = 191.6301;
A_frac_calc = ((C11_max-1) - (C11_min-1))/(2*(L-1)) %*(L-1)/2
%}
%%
function [Q_value] = calc_Q_value(L, fN, fnn, pw, N_inf, f_MF)
    % N_inf = 1: N->infty limit
    % mean-field contribution
    if ~N_inf
        f_MF = fN'-6*fnn;
    end
    
    % calculate intersection points and define integrals
    syms Con % for integrals
    % ---self interaction---
    % intersection points 
    Cmax_self = (L - 1 - 4.*fnn - f_MF.*(1-pw))./(2.*fnn+f_MF.*pw);
    Cmin_self = (L - 1 - 6.*fnn - f_MF.*(1-pw))./(f_MF.*pw);
    
    % Integrals
    I_self = cell(2, 1);
    I_self{1} = @(Con) 2.*fnn.*(Con - 1);
    I_self{2} = @(Con) (L - 1 - 6.*fnn - f_MF.*(Con.*pw + 1-pw));
    I_self{3} = @(Con) 0;

    % ---mutual interaction---
    % intersection points 
    Cmax_mutual = (L - 0 - 2.*fnn - f_MF.*(1-pw))./(1 + 4.*fnn+f_MF.*pw);
    Cmin_mutual = (L - 1 - 6.*fnn - f_MF.*(1-pw))./(2*fnn + f_MF.*pw);

    % Integrals
    I_mutual = cell(2, 1);
    I_mutual{1} = @(Con) (1+2.*fnn).*(Con - 1);
    I_mutual{2} = @(Con) (L - 1 - 2.*Con.*fnn - 4.*fnn - f_MF.*(Con.*pw + 1-pw));
    I_mutual{3} = @(Con) 0;

    %% Interactions 1 <- 1 and 2 <- 1
    % determine integral bounds
    [Con_bounds, K_bounds] = get_all_bounds(Cmax_self(1), Cmin_self(1), Cmax_mutual(1), Cmin_mutual(1), L);
    Con_bound_vals = [1 Con_bounds];
    
    %disp('Con_bound_vals = ');
    %disp(Con_bound_vals); %-------------
    % disp([Ca_self(1), Cb_self(1), Ca_mutual(1), Cb_mutual(1)]);
    
    % calculate integral
    mol = 1; % molecule 1
    I_sum = 0; % integral sum value
    syms x
    for i=1:numel(Con_bounds)
        % calculate integral over bounds
        integrand_term1 = I_self{K_bounds(1, i)}(x);
        integrand_term2 = I_mutual{K_bounds(2, i)}(x);
        if ~(integrand_term1(mol)==0 && integrand_term2(mol)==0)
            this_I = double( int(integrand_term1(mol)*integrand_term2(mol), Con_bound_vals(i), Con_bound_vals(i+1)) );
            I_sum = I_sum + this_I;
        end
    end
    frac_mol_1 = I_sum/(L-1)^3;
    %% Interaction 1 <- 2 (new method)
    Ca = Cmax_mutual(2);
    Cb = Cmin_mutual(2);
    mol = 2;
    
    %disp('Ca Cb = ');
    %disp([Ca Cb]); %-------------
    
    if Ca<L && Cb<L
        %disp('case A');
        integrand_1 = I_mutual{1}(x);
        I1 =  double( int(integrand_1(mol), 1, Ca) );
        integrand_2 = I_mutual{2}(x);
        I2 =  double( int(integrand_2(mol), Ca, Cb) );
        frac_mol_2 = (I1 + I2)/(L-1)^2;
    elseif Ca<L && Cb>L
        %disp('case B');
        integrand_1 = I_mutual{1}(x);
        I1 =  double( int(integrand_1(mol), 1, Ca) );
        integrand_2 = I_mutual{2}(x);
        I2 =  double( int(integrand_2(mol), Ca, L) );
        frac_mol_2 = (I1 + I2)/(L-1)^2;
    else
        %disp('case C');
        integrand_1 = I_mutual{1}(x);
        I1 = double( int(integrand_1(mol), 1, L) );
        frac_mol_2 = I1/(L-1)^2;
    end

    % Final result
    %disp(frac_mol_1)
    %disp(frac_mol_2)
    Q_value = frac_mol_1*frac_mol_2;
end

function [fN, fnn, fnnn] = calc_fN(a0, rcell, gz, dist, lambda)
    % calculate fN
    Rcell = a0*rcell;
    idx = gz + round(gz/2); % pick cell not at corner of grid
    dist_vec = a0*dist(idx,:);
    r = dist_vec(dist_vec>0);
    fN = zeros(2,1);
    fN(1) = sum(sinh(Rcell./lambda(1))*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
    fN(2) = sum(sinh(Rcell./lambda(2))*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));

    % calculate fnn
    fnn = zeros(1, 2);
    rnn = a0;
    fnn(1) = sinh(Rcell./lambda(1))*exp((Rcell-rnn)./lambda(1)).*(lambda(1)./rnn) ; % calculate signaling strength
    fnn(2) = sinh(Rcell./lambda(2))*exp((Rcell-rnn)./lambda(2)).*(lambda(2)./rnn) ; % calculate signaling strength

    fnnn = zeros(2,2);
    rnnn = [sqrt(3); 2].*a0;
    fnnn(1,:) = sinh(Rcell./lambda(1))*exp((Rcell-rnnn)./lambda(1)).*(lambda(1)./rnnn);
    fnnn(2,:) = sinh(Rcell./lambda(2))*exp((Rcell-rnnn)./lambda(2)).*(lambda(2)./rnnn);
end

function [Con_bounds, K_bounds] = get_all_bounds(C1a, C1b, C2a, C2b, L)
    % C1a from Kmax = L for molecule 1
    % C1b from Kmin = L for molecule 1 
    % C2a from Kmax = L for molecule 2 
    % C2b from Kmin = L for molecule 2 
    
    % returns integral_bounds: 2 x n matrix, n<=5
    % specifies boundaries of each of the n regions for both molecules (see get_bounds)
    Con_bounds = sort([C1a C1b C2a C2b L]);
    Con_bounds = Con_bounds(Con_bounds <= L);

    K_bounds = zeros( 2, numel(Con_bounds));
    for i=1:numel(Con_bounds)
        upper_bound = Con_bounds(i);
        K_bounds(1, i) = get_bounds(C1a, C1b, upper_bound);
        K_bounds(2, i) = get_bounds(C2a, C2b, upper_bound);
    end
end
%squeeze( integral_bounds(1, :, :) )
%squeeze( integral_bounds(2, :, :) )

function bounds = get_bounds(Ca, Cb, b1)
    % 3 levels: convert to numbers    
    if Ca >= b1 && Cb >= b1
        bounds = 1; %case 1: [Kmin Kmax];
    elseif Ca <= b1 && Cb >= b1
        bounds = 2; %case 2: [Kmin L];
    else
        bounds = 3; %case 3:[L L];
    end
end

function h = plot_boundaries_K_Con(fN, fnn, pw, L, showlines)
    Con1_all = linspace(1, L, 100);
    Con2_all = linspace(1, L, 100);
    YMF1_all = (fN(1) - 6*fnn(1))*(Con1_all*pw(1) + (1-pw(1)));
    YMF2_all = (fN(2) - 6*fnn(2))*(Con2_all*pw(2) + (1-pw(2)));

    % bounds
    K11_all_upb = 1 + 2*Con1_all*fnn(1) + 4*fnn(1) + YMF1_all;
    K11_all_lwb = 1 + 6*fnn(1) + YMF1_all;
    K12_all_upb = Con2_all + (4*Con2_all + 2)*fnn(2) + YMF2_all;
    K12_all_lwb = 1 + (2*Con2_all + 4)*fnn(2) + YMF2_all;
    K21_all_upb = Con1_all + 4*Con1_all*fnn(1) + 2*fnn(1) + 4*fnn(1) + YMF1_all;
    K21_all_lwb = 1 + 2*Con1_all*fnn(1) + 4*fnn(1) + YMF1_all;

    % intersection points 
    C1_self = (L - 1 - 4.*fnn - (fN' - 6*fnn).*(1-pw))./(2.*fnn+(fN'-6.*fnn).*pw);
    C2_self = (L - 1 - 6.*fnn - (fN' - 6*fnn).*(1-pw))./((fN'-6.*fnn).*pw);

    C1_mutual = (L - 0 - 2.*fnn - (fN' - 6*fnn).*(1-pw))./(1 + 4.*fnn+(fN'-6.*fnn).*pw);
    C2_mutual = (L - 1 - 6.*fnn - (fN' - 6*fnn).*(1-pw))./(2*fnn + (fN'-6.*fnn).*pw);

    h=figure;
    subplot(2,2,1);
    hold on
    plot(K11_all_upb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
    plot(K11_all_lwb, Con2_all, '--', 'Color', 'r', 'LineWidth', 2 )
    if showlines
        plot([1 L], [C1_self(1) C1_self(1)], '--', 'Color', 'b');
        plot([1 L], [C2_self(1) C2_self(1)], '--', 'Color', 'b');
    end
    xlim([0 L]);
    ylim([0 L]);
    title('1 \leftarrow 1');
    xlabel('K^{(11)}');
    ylabel('C_{ON}^{(1)}');
    set(gca, 'FontSize', 20);

    subplot(2,2,2);
    hold on
    plot(K12_all_upb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
    plot(K12_all_lwb, Con2_all, '--', 'Color', 'r', 'LineWidth', 2 )
    if showlines
        plot([1 L], [C1_mutual(2) C1_mutual(2)], '--', 'Color', 'b');
        plot([1 L], [C2_mutual(2) C2_mutual(2)], '--', 'Color', 'b');
    end
    xlim([0 L]);
    ylim([0 L]);
    title('1 \leftarrow 2');
    xlabel('K^{(12)}');
    ylabel('C_{ON}^{(2)}');
    set(gca, 'FontSize', 20);

    subplot(2,2,3);
    hold on
    plot(K21_all_upb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
    plot(K21_all_lwb, Con2_all, '--', 'Color', 'r', 'LineWidth', 2 )
    if showlines
        plot([1 L], [C1_mutual(1) C1_mutual(1)], '--', 'Color', 'b');
        plot([1 L], [C2_mutual(1) C2_mutual(1)], '--', 'Color', 'b');
    end
    xlabel('K^{(21)}');
    ylabel('C_{ON}^{(1)}');
    xlim([0 L]);
    ylim([0 L]);
    title('2 \leftarrow 1');
    set(gca, 'FontSize', 20);
end

function [A, C1, C2] = calc_area(L, fN, fnn, pw, interaction, molecule, N_inf, f_MF)
    % interaction: 'self' -> self (1); 'mutual' -> mutual (2)
    % molecule: 1 or 2    
    % N_inf = 1: N->infty limit
    
    % convert string input to numeric
    int = find( cellfun(@(x) ~isempty(x), regexp( interaction, {'self', 'mutual'}, 'match' )), 1);
    
    % specify molecule 
    fN = fN(molecule);
    fnn = fnn(molecule);
    pw = pw(molecule);
    
    % mean-field approximation?
    if ~N_inf % if not N-> infty
        f_MF = (fN'-6.*fnn); % restore original expression
    end
    
    % calculate intersection points and define integrals
    syms x x1 x2 % for integrals
    switch int
        case 1 % self
            % intersection points 
            C1 = (L - 1 - 4.*fnn - f_MF.*(1-pw))./(2.*fnn+f_MF.*pw);
            C2 = (L - 1 - 6.*fnn - f_MF.*(1-pw))./(f_MF.*pw);

            % Integrals
            I1 = @(x) (x.^2 - 2*x - 1).*fnn;
            I2 = @(x1, x2) L*(x2-x1) - f_MF.*pw/2.*(x2.^2-x1.^2) -...
                (1 + 6*fnn + f_MF.*(1-pw)).*(x2-x1);
        case 2 % mutual
            % intersection points 
            C1 = (L - 0 - 2.*fnn - f_MF.*(1-pw))./(1 + 4.*fnn+f_MF.*pw);
            C2 = (L - 1 - 6.*fnn - f_MF.*(1-pw))./(2*fnn + f_MF.*pw);

            % Integrals
            I1 = @(x) (x.^2 - 2*x + 1)/2.*(1+2*fnn); % integral 1
            I2 = @(x1, x2) L*(x2-x1) - ( (fnn + f_MF.*pw/2).*(x2.^2-x1.^2) -...
                (1 + 4*fnn + f_MF.*(1-pw)).*(x2-x1) ) ; % integral 2 -> check calculation
    end

    % Boundary cases
    if C1<L && C2<L
        % case A
        %disp('case A');
        A = I1(C1) + I2(C1, C2);
    elseif C1<L && C2>L
        % case B   
        %disp('case B');
        A = I1(C1) + I2(C1, L);
    else
        % case C
        %disp('case C');
        A = I1(L);
    end
    A = A/(L-1)^2;
end