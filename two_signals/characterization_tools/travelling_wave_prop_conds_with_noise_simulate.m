% Calculate the TW propagation conditions with noise
clear all
close all
set(0,'defaulttextinterpreter', 'tex')

%% Parameters
% Manual input
network = 36;
%M_int = [1 -1; 1 0]; % network 15
%M_int = [1 1; -1 0]; % network 19
%M_int = [1 -1; 1 1]; % network 33
%M_int = [-1 -1; 1 1]; % network 34
M_int = [-1 1; -1 1]; % network 36

gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];
hill = Inf;
Coff = [1 1];

% K, Con (crucial)
Con = [500 500];
K = [950 400; 200 80];

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% calculate fN
[fN, fnn, fnnn] = calc_fN(a0, rcell, gz, dist, lambda);

% Simulation parameters
alpha_all_sim = 10.^[-2:0.1:0];
num_sim = 100;

run = 7;
IO_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_robustness_extensions\vs_noise\run_7_network_36';

%% Check boundaries in K, Con space (only network 15)
%{
pw = [1/8 1/8];
L = 1000;
showlines = 0;
h = plot_K_Con_with_boundaries(fN, fnn, pw, L, showlines, K, Con);
%}
%% ------------------------------------------------------------------------
% Part II: Simulations
% Check that wave can indeed propagate for the zero noise condition
noise = 0;

% Load and plot initial state
show_fig = 0;
[~, cells_ini] = plot_ini_state(network, gz, rcell, a0, show_fig);

% Run simulation
display_fig = 1;
sim_ID = 'two_signal_mult';
mcsteps = 0;
InitiateI = 0;
p0 = 0; 
I0 = 0;
tmax = 10^4;
save_folder = '';
fname_str_template = strrep(sprintf('%s_N%d_ini_state_TW_noise_%.3f',...
    sim_ID, N, noise), '.', 'p');

[cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks(...
    N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
    dist, pos, sim_ID, mcsteps, InitiateI, p0, I0, cells_ini,...
    tmax, save_folder, fname_str_template, display_fig);

%% Batch simulations
%
% determine # sims to do    
sim_to_do = calc_num_sim_to_do(IO_folder, alpha_all_sim, sim_ID, N, num_sim); % ! check filename pattern

% Do simulations; loop over noise values
for idx_param_loop=1:numel(alpha_all_sim)
    noise = alpha_all_sim(idx_param_loop);
    
    sim_to_do_temp = sim_to_do(idx_param_loop);
    for trial=1:sim_to_do_temp
        fprintf('Noise %.3f, trial %d/%d \n', noise, trial, sim_to_do_temp);
        
        display_fig = 0;
        fname_str_template = strrep(sprintf('%s_N%d_ini_state_TW_params_noise_%.3f',...
        	sim_ID, N, noise), '.', 'p');
        % previously, time_evolution_save_func_efficient_checks
        [cells_hist, period, t_onset] = check_TW_propagation_save(...
            N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
            dist, pos, sim_ID, mcsteps, InitiateI, p0, I0, cells_ini,...
            tmax, IO_folder, fname_str_template, display_fig); 
    end 
end
%
%% Analyze simulation results
pattern = strrep(sprintf('%s_N%d_ini_state_TW_params_noise_%s_t_out_%s_period_%s',...
	sim_ID, N, '(\d)p(\d+)', '(\d+)', '(\d+|Inf)' ), '.', 'p');        

% get all filenames
listing = dir(IO_folder);
num_files = numel(listing)-2;
names = {};
filecount = 0;
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        match = regexp(name, pattern, 'match');
        %disp(match);
        if ~isempty(match)
            filecount = filecount + 1;
            names{end+1} = name;
        end
    end
end

%%
% Load data and store results
TW_persistence = zeros(numel(alpha_all_sim), num_sim);
TW_persistence_2 = zeros(numel(alpha_all_sim), num_sim);
sim_count = zeros( numel(alpha_all_sim), 1 );

for idx=1:numel(names)
    load( fullfile(IO_folder, names{idx}), 'cells_hist', 't_out', 'period',...
        'save_consts_struct');
    disp(names{idx});
    idx_noise = find(round(alpha_all_sim, 3) == round(save_consts_struct.noise, 3));
    sim_count(idx_noise) = sim_count(idx_noise) + 1;
    
    % determine whether TW persisted or not
    digits = 3;
    [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
        period, t_out, dist, digits);
    TW_persistence(idx_noise, sim_count(idx_noise)) = trav_wave;
    TW_persistence_2(idx_noise, sim_count(idx_noise)) = trav_wave_2;   
    disp(trav_wave_2);
end
%
data_sim = sum(TW_persistence_2, 2)/num_sim;

%% Save analyzed results
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_robustness_extensions\vs_noise';
fname_str = sprintf('analyzed_TW_robustness_vs_noise_run_%d', run);

save_vars = {N, a0, K, Con, Coff, M_int, hill,...
    rcell, Rcell, lambda,...
    sim_ID, tmax, mcsteps};
save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill',...
    'rcell', 'Rcell', 'lambda', ...
    'sim_ID', 'tmax', 'mcsteps'};
save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);

save( fullfile(save_folder, fname_str), 'alpha_all_sim', 'TW_persistence',...
    'TW_persistence_2', 'save_consts_struct' );
%--------------------------------------------------------------------------
%}

%% Load analyzed results
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_robustness_extensions\vs_noise';
fname_str = sprintf('analyzed_TW_robustness_vs_noise_run_%d', run);
load( fullfile(save_folder, fname_str), 'alpha_all_sim', 'TW_persistence',...
    'TW_persistence_2', 'save_consts_struct' );
data_sim = sum(TW_persistence_2, 2)/num_sim;

%% Plot analyzed results
%alpha_all_sim = 10.^(-2:0.5:0);
h=figure;
semilogx(alpha_all_sim, data_sim, 'bo-');
xlabel('Noise strength \alpha');
ylabel('Fraction TW persisted');
title('TW persistence');
set(gca, 'FontSize', 24);

% Save plot
qsave = 1;
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_robustness_extensions\vs_noise';
fname_str = sprintf('analyzed_TW_robustness_vs_noise_run_%d_plot_sim', run);
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Calculate theoretical results
alpha_all_calc = 10.^[-2:0.01:0];
%--------------------------------------------------------------------------
% state permutations of wave
switch network
    case 15
        %states_perm = [2 4 3 1]; % network 15 reversed
        states_perm = [3 4 2 1]; % network 15 
    case 19
        states_perm = [4 3 1 2]; % network 19
    case 33 % ! adjust manually wave form !
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

P_propagation_all = zeros( numel(alpha_all_calc), 1 );
for i=1:numel(alpha_all_calc)
    alpha = alpha_all_calc(i);
    P_propagation = calc_P_prop(Y_all, target_states_all, M_int, K, alpha, gz);
    P_propagation_all(i) = P_propagation;
end

% plot theoretical results only (check)
%h=figure;
%semilogx(alpha_all_calc, P_propagation_all, 'k-');

%% Plot analyzed results together with theoretical results
% Plot TW survival (theory + simulations)
h=figure;
set(gca, 'XScale', 'log');
hold on
semilogx(alpha_all_calc, P_propagation_all, 'k-');
%semilogx(alpha_all_calc, P_survival_all, 'k--');
scatter(alpha_all_sim, data_sim, 64, 'ko');
xlabel('Noise strength \alpha');
ylabel('Survival probability');
set(gca, 'FontSize', 24);
legend({'Theory', 'Simulations'}, 'Location', 'ne');
%legend({'Theory (full)', 'Theory (approx)', 'Simulations'}, 'Location', 'ne');

% Save plot
qsave = 1;
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_robustness_extensions\vs_noise';
fname_str = sprintf('analyzed_TW_robustness_vs_noise_run_%d_plot_theory_and_sim', run);
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% --- "Simple way of calculating P_survival (not accurate)" ---
% Calculate P_survival
%alpha_all_calc = 10.^[-2:0.01:0]; %[0.01 0.05 0.1 0.5 1 5 10];
%{
% Calculate bounds of K^{i,j}
pw = [1/8 1/8]; % mean-field contribution
[dK_upper_bounds, dK_lower_bounds] = calc_K_prop_bounds(K, Con, fN, fnn, pw);
P_survival_all = zeros(size(alpha_all_calc));
for i=1:numel(alpha_all_calc)
    alpha = alpha_all_calc(i);
    P_aux = normcdf(dK_upper_bounds, 0, alpha.*K)-normcdf(dK_lower_bounds, 0, alpha.*K);
    P_survival_all(i) = prod(P_aux(abs(M_int)>0));
end

% Plot TW survival (theory + simulations)
h=figure;
semilogx(alpha_all_calc, P_survival_all, '-');
hold on
semilogx(alpha_all_sim, data_sim, 'o');
xlabel('Noise strength \alpha');
ylabel('Survival probability');
set(gca, 'FontSize', 24);
%% ------------------------ Troubleshooting -------------------------------

alpha = 0.05;
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------   
P_target_all = zeros(6, 2); % Probability that with noise, the system reaches the right target
% idx = 1; 
for idx = 3 %1:6 % index of cell (E_F, F, M, etc.)
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

%% Sim data
% Load and plot initial state
show_fig = 0;
[~, cells_ini] = plot_ini_state(network, gz, rcell, a0, show_fig);

idx_M = find(cells_ini*[2; 1] == 3);
idx_B = find(cells_ini*[2; 1] == 1);
idx_F = find(cells_ini*[2; 1] == 2);

% Calculate Y_sensed
[cells_out, changed, Y] = ...
    update_cells_two_signals_multiply_finite_Hill(cells_ini, dist, M_int, a0,...
    Rcell, Con, Coff, K, lambda, hill, noise);
Y(idx_M, :) % Y_sensed of M cells
Y(idx_B, :) % Y_sensed of B cells
%Y(idx_F, :) % Y_sensed of B cells
%}

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

function sim_to_do = calc_num_sim_to_do(folder, noise_all, sim_ID, N, num_sim)
    if exist(folder, 'dir') ~= 7
        warning('Folder does not exist! ');
        mkdir(folder);
        fprintf('Made new folder %s \n', folder);
    end
    
    for idx_loop=1:numel(noise_all)
        noise = noise_all(idx_loop);
        %(!!!)  % Filename pattern (!!!)
        pattern = strrep(sprintf('%s_N%d_ini_state_TW_params_noise_%.3f_t_out_%s_period_%s',...
            sim_ID, N, noise, '(\d+)', '(\d+|Inf)' ), '.', 'p');

        listing = dir(folder);
        num_files = numel(listing)-2;
        names = {};
        filecount = 0;
        for i = 1:num_files
            filename = listing(i+2).name;
            % remove extension and do not include txt files
            [~,name,ext] = fileparts(filename);
            if strcmp(ext, '.mat')
                match = regexp(name, pattern, 'match');
                %disp(match);
                if ~isempty(match)
                    filecount = filecount + 1;
                    names{end+1} = name;
                end
            end
        end
        fprintf('N=%d, noise = %.3f, sim to do: %d \n',...
            N, noise, num_sim-filecount);
        sim_to_do(idx_loop) = num_sim-filecount;
    end
    fprintf('Total number of simulations to do: %d \n', sum(sim_to_do(:)) );
end

% --------------- Theoretical ---------------------------------------------
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

%-------------- Discarded (inaccurate/wrong)----------------------------

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

%----------------------------- theory: plot region ------------------------

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
