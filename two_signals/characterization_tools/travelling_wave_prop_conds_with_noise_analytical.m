% Calculate the TW propagation conditions with noise
clear all
close all
set(0,'defaulttextinterpreter', 'tex')

%% Parameters
% Manual input
network = 15;
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];
hill = Inf;
Coff = [1 1];

% K, Con (crucial)
M_int = [1 -1; 1 0];
Con = [500 500];
K = [50 300; 300 0];

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% calculate fN
[fN, fnn, fnnn] = calc_fN(a0, rcell, gz, dist, lambda);

% mean-field contribution
pw = [1/8 1/8];
L = 1000;
Con1_all = linspace(1, L, 100);
Con2_all = linspace(1, L, 100);
YMF1_all = (fN(1) - 6*fnn(1))*(Con1_all*pw(1) + (1-pw(1)));
YMF2_all = (fN(2) - 6*fnn(2))*(Con2_all*pw(2) + (1-pw(2)));
%% ------------------------------------------------------------------------
% Part I: Theory
% Plot parameter set together with propagation conditions
showlines = 0;
h = plot_K_Con_with_boundaries(fN, fnn, pw, L, showlines, K, Con);

% Calculate allowed bounds for K fluctuations
% Calculate bounds of K^{i,j}
[dK_upper_bounds, dK_lower_bounds] = calc_K_prop_bounds(K, Con, fN, fnn, pw);

% Calculate survival probability 
% Single instance - based on final inequalities
%{
alpha = 1; 
P_aux = normcdf(dK_upper_bounds, 0, alpha.*K)-normcdf(dK_lower_bounds, 0, alpha.*K);
P_survival = prod(P_aux(abs(M_int)>0));
fprintf('Survival probability = %.5f \n', P_survival);

% Versus alpha
alpha_all = 10.^[-2:0.01:1]; %[0.01 0.05 0.1 0.5 1 5 10];
P_survival_all = zeros(size(alpha_all));
for i=1:numel(alpha_all)
    alpha = alpha_all(i);
    P_aux = normcdf(dK_upper_bounds, 0, alpha.*K)-normcdf(dK_lower_bounds, 0, alpha.*K);
    P_survival_all(i) = prod(P_aux(abs(M_int)>0));
end
%}

%% Based on explicit conditions
alpha = 0;
%--------------------------------------------------------------------------
% state permutations of wave
switch network
    case 15
        %states_perm = [2 4 3 1]; % network 15 reversed
        states_perm = [3 4 2 1]; % network 15 
    case 19
        states_perm = [4 3 1 2]; % network 19
    case 33
        states_perm = [3 4 2 1]; % network 33/15
        %states_perm = [4 2 1 3]; % network 33/34
    case 34
        states_perm = [4 2 1 3]; % network 33/34
    case 36
        states_perm = [2 4 3 1]; % network 36
end
default_states = [0 0; 0 1; 1 0; 1 1];
states = default_states(states_perm, :); % F, M, B, E
%--------------------------------------------------------------------------
% Calculate sensed concentrations
Y_all = calc_Y_all(gz, Con, fN, fnn, states_perm);
%--------------------------------------------------------------------------
% determine output states
% (1) E_F->F, (2) F->M, (3) M->B, (4) B->E, (5) E_B->E, (6) E->E
targets_all = [1 2 3 4 4 4]; % recall F, M, B, E
target_states_all = states(targets_all,:);
%% Calculate propagation probability (for 1 time step)
% Noise conditions
alpha_all = 10.^[-2:0.1:0];
P_propagation_all = zeros( numel(alpha_all), 1 );
for i=1:numel(alpha_all)
    alpha = alpha_all(i);
    P_propagation = calc_P_prop(Y_all, target_states_all, M_int, K, alpha, gz);
    P_propagation_all(i) = P_propagation;
end

% Plot figure
h = figure;
semilogx(alpha_all, P_propagation_all, 'o-');
xlabel('Noise strength \alpha');
ylabel('Survival probability');
title('Survival after one time step');
set(gca, 'FontSize', 24);

%% Plot survival probability after x time steps
ts = [1 15]; %10 100 10^5];
figure;
for i=1:numel(ts)
    semilogx(alpha_all, P_propagation_all.^ts(i), 'x-');
    hold on
end
xlabel('Noise strength \alpha');
ylabel('Survival probability');
set(gca, 'FontSize', 24);
legend(sprintfc('%d', ts));


%% Functions
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


function h = plot_K_Con_with_boundaries(fN, fnn, pw, L, showlines, K, Con)
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

    h=figure; % 1<-1
    subplot(2,2,1);
    hold on
    plot(K11_all_upb, Con1_all, '-', 'Color', 'r', 'LineWidth', 2 )
    plot(K11_all_lwb, Con1_all, '--', 'Color', 'r', 'LineWidth', 2 )
    if showlines
        plot([1 L], [C1_self(1) C1_self(1)], '--', 'Color', 'b');
        plot([1 L], [C2_self(1) C2_self(1)], '--', 'Color', 'b');
    end
    % parameter set
    plot(K(1,1), Con(1), 'bx')
    xlim([0 L]);
    ylim([0 L]);
    title('1 \leftarrow 1');
    xlabel('K^{(11)}');
    ylabel('C_{ON}^{(1)}');
    set(gca, 'FontSize', 20);

    subplot(2,2,2); %1<-2
    hold on
    plot(K12_all_upb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
    plot(K12_all_lwb, Con2_all, '--', 'Color', 'r', 'LineWidth', 2 )
    if showlines
        plot([1 L], [C1_mutual(2) C1_mutual(2)], '--', 'Color', 'b');
        plot([1 L], [C2_mutual(2) C2_mutual(2)], '--', 'Color', 'b');
    end
    % parameter set
    plot(K(1,2), Con(2), 'bx')
    xlim([0 L]);
    ylim([0 L]);
    title('1 \leftarrow 2');
    xlabel('K^{(12)}');
    ylabel('C_{ON}^{(2)}');
    set(gca, 'FontSize', 20);

    subplot(2,2,3); %2<-1
    hold on
    plot(K21_all_upb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
    plot(K21_all_lwb, Con2_all, '--', 'Color', 'r', 'LineWidth', 2 )
    if showlines
        plot([1 L], [C1_mutual(1) C1_mutual(1)], '--', 'Color', 'b');
        plot([1 L], [C2_mutual(1) C2_mutual(1)], '--', 'Color', 'b');
    end
    % parameter set
    plot(K(2,1), Con(1), 'bx')
    xlabel('K^{(21)}');
    ylabel('C_{ON}^{(1)}');
    xlim([0 L]);
    ylim([0 L]);
    title('2 \leftarrow 1');
    set(gca, 'FontSize', 20);
end

function [dK_upper_bounds, dK_lower_bounds] = calc_K_prop_bounds(K, Con, fN, fnn, pw)
    dK_upper_bounds = zeros(2);
    dK_lower_bounds = zeros(2);

    YMF1 = (fN(1) - 6*fnn(1))*(Con(1)*pw(1) + (1-pw(1)));
    YMF2 = (fN(2) - 6*fnn(2))*(Con(2)*pw(2) + (1-pw(2)));

    K11_upb = 1 + 2*Con(1)*fnn(1) + 4*fnn(1) + YMF1;
    K11_lwb = 1 + 6*fnn(1) + YMF1;
    K12_upb = Con(2) + (4*Con(2) + 2)*fnn(2) + YMF2;
    K12_lwb = 1 + (2*Con(2) + 4)*fnn(2) + YMF2;
    K21_upb = Con(1) + 4*Con(1)*fnn(1) + 2*fnn(1) + 4*fnn(1) + YMF1;
    K21_lwb = 1 + 2*Con(1)*fnn(1) + 4*fnn(1) + YMF1;

    dK_upper_bounds(1,1) = K11_upb - K(1,1);
    dK_upper_bounds(1,2) = K12_upb - K(1,2);
    dK_upper_bounds(2,1) = K21_upb - K(2,1);
    dK_lower_bounds(1,1) = K11_lwb - K(1,1);
    dK_lower_bounds(1,2) = K12_lwb - K(1,2);
    dK_lower_bounds(2,1) = K21_lwb - K(2,1);
end

function [h, cells_ini] = plot_ini_state(network, gz, rcell, a0, show_fig)
    N = gz^2;
    % get network_idx
    networks_all = [15 19 33 33 34 36]; 
    appendix = ''; % Note special rule for 33a, 33b
    network_idx = find(network==networks_all, 1);
    if strcmp(appendix, 'b')
        network_idx = 4; % special case: 33b
    end

    signal_count = 2;
    folder = 'H:\My Documents\Multicellular automaton\app\data\system_states';
    fname = fullfile(folder, 'trav_wave_single_vertical_central_position');
    [~, cells_ini, ~] = manual_input_state(signal_count, folder, N, fname);

    cell_states_all = cell(6, 1);
    cell_states_all{1} = [1 0; 1 1; 0 1; 0 0]; % F, M, B, E
    cell_states_all{2} = [1 1; 1 0; 0 0; 0 1]; % F, M, B, E
    cell_states_all{3} = cell_states_all{1}; % F, M, B, E
    cell_states_all{4} = [1 1; 0 1; 0 0; 1 0]; % F, M, B, E
    cell_states_all{5} = cell_states_all{4}; % F, M, B, E
    cell_states_all{6} = [0 1; 1 1; 1 0; 0 0]; % F, M, B, E
    cell_states = cell_states_all{network_idx};

    cells_idx00_E = find( cells_ini*[1; 2]==0 );
    cells_idx01_F = find( cells_ini*[1; 2]==1 );
    cells_idx10_B = find( cells_ini*[1; 2]==2 );
    cells_idx11_M = find( cells_ini*[1; 2]==3 );

    cells_ini(cells_idx01_F, :) = repmat(cell_states(1,:), gz, 1);
    cells_ini(cells_idx11_M, :) = repmat(cell_states(2,:), gz, 1);
    cells_ini(cells_idx10_B, :) = repmat(cell_states(3,:), gz, 1);
    cells_ini(cells_idx00_E, :) = repmat(cell_states(4,:), N-3*gz, 1);

    % show lattice
    if show_fig
        nodisplay = 1; 
        mcsteps = 0;
        [pos_ini, dist_ini] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);

        % check initial state
        h = figure;
        plot_handle = reset_cell_figure(h, pos_ini, rcell);
        t=0; disp_mol = 12; showI = 0; 
        update_figure_periodic_scatter(plot_handle, cells_ini, t, disp_mol, showI, a0, dist_ini);
    else
        h = [];
    end
end

function Y_all = calc_Y_all(gz, Con, fN, fnn, states_perm)
    % Calculate Y_all
    default_states = [0 0; 0 1; 1 0; 1 1];
    
    % calculate Y(alpha)
    % neighbour data
    n_nei = [2	0	0	4; % EF
        2	2	0	2; % F
        2	2	2	0; % M
        0	2	2	2; % B
        0	0	2	4; % EB
        0	0	0	6]; % E
    states = default_states(states_perm, :); % F, M, B, E
    types = [states(4,:); states(1,:); states(2,:); states(3,:); states(4,:); states(4,:)];
    Y_self = (Con-1).*types + 1;

    % calculate Y_nei
    Y_nei = zeros(size(n_nei, 1), 2);
    for i=1:6
        Y_nei(i,1) = fnn(1)*n_nei(i,:)*((Con(1)-1)*states(:,1)+1);
        Y_nei(i,2) = fnn(2)*n_nei(i,:)*((Con(2)-1)*states(:,2)+1);
    end

    % estimate p
    tmp = 1; %num_waves*bandwidth;
    tmp2 = [tmp tmp tmp gz-3*tmp];
    p = (tmp2*states)/gz;

    % calculate Y_mf (mean-field contribution)
    z = 6; % coordination number
    Y_mf = (fN' - z*fnn).*( Con.*p + (1-p) );

    Y_all = Y_self + Y_nei + Y_mf;
end

function [P_propagation, P_target_all] = calc_P_prop(Y_all, target_states_all, M_int, K, alpha, gz)
    P_target_all = zeros(6, 2); % Probability that with noise, the system reaches the right target
    % idx = 1; 
    for idx = 1:6 % index of cell (E_F, F, M, etc.)
        % set sensed conc. and target for this wave state
        Y_idx = Y_all(idx,:);
        target = target_states_all(idx, :);
        % probability( g^{ij}_alpha = 1 ), i.e. the interaction is unrepressed
        P_unrepressed = abs(M_int).*((1+M_int)/2.*normcdf(Y_idx - K, 0, alpha.*K) + (1-M_int)/2.*(1-normcdf(Y_idx - K, 0, alpha.*K)) ); 
        P_unrepressed(abs(M_int)==0) = 1; % absent interactions: 
        % probability( output = target )
        P_target = target.*prod(P_unrepressed, 2)' + (1-target).*(1 - prod(P_unrepressed, 2))';
        P_target_all(idx, :) = P_target;
    end
    N_wave_state = [gz; gz; gz; gz; gz; (gz-5)*gz]; % number of cells with given state
    P_propagation = prod(prod(P_target_all, 2).^N_wave_state);
end