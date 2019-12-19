% Analyze parameter sets for which TWs can propagate in more detail
% Plot K-Con plots for each interaction, compare 
% Case study: network 15
% Wrong calculation due to product of areas
% Plot code still useful
clear all
close all
set(0, 'defaulttextinterpreter', 'tex')

%% Parameters
% Number of parameter sets to do
n_pset = 10^5;

% Manual input
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
%Rcell = rcell*a0;
lambda = [1 1.2];
Coff = [1 1];

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% specify wave type and characteristics
wave_types_str = {'straight', 'inward bend', 'outward bend'};
wave_type = 1;
num_waves = 1; % number of waves
bandwidth = 1; % width of band of cells of one type

% Data folder
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
subfolder = 'run2';

% save figure folder
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\robustness_analytical';

%% Load analyzed data on TW propagation (analytical)
%
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_%s',...
    num_waves, wave_type, subfolder);
load( fullfile(folder, fname_str) );
%}

%% display found waves
%
[x_found, y_found] = find(wave_possible);
t=table(P(x_found, :), y_found, 'VariableNames', {'Wave_type', 'Network'});
t2=table(x_found, y_found, 'VariableNames', {'wave_idx', 'Network'});

disp('Cell states F, M, B, E');
disp(t);
disp(t2);

% Get a list of all interaction matrices for the found networks
[M_int_found] = get_found_M_int(y_found);
%}

% Network 15 settings
idx_loop = 1;
wave_idx = x_found(idx_loop);
network = y_found(idx_loop);
fprintf('wave_idx %d, network %d \n', wave_idx, network);
states_perm = P(wave_idx, :);

% calculate fN, fnn
%[fN, fnn, ~] = calc_fN(a0, rcell, gz, dist, lambda);
% pw = [2/gz 2/gz];
%% Load and process detailed data per network
%{
%for idx_loop=1 %1:numel(x_found)
    idx_loop = 1;
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    states_perm = P(wave_idx, :);
    %
    %--------------------------------------------------------------------------
    % Get K, Con parameters for waves (analytical)
    % load data
    %folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
    %subfolder = 'run3_vary_a0_lambda12';
    fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
            num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    load( fullfile(folder, subfolder, fname_str) );
    
    Con_all_network = squeeze(Con_all(network, :, :));
    K_all_network = squeeze(K_all(network, :, :, :));
    trav_wave_conds_met_theory = squeeze(trav_wave_conds_met(network, :, :))';
    
    Con_wave_theory = squeeze(Con_all(network, trav_wave_conds_met_theory, :));
    K_wave_theory = squeeze(K_all(network, trav_wave_conds_met_theory, :, :));
    num_wave_sims = size(Con_wave_theory, 1);
    
    clear('Con_all', 'K_all', 'trav_wave_conds_met')
    %M_int = M_int_found{idx_loop}; 
    %%
    % Get K, Con parameters for waves (simulations)
    % load data
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2_stability_sim';
    %if remote
    %    load_folder = strrep(load_folder, 'N:\', 'W:\staff-bulk\');
    %end
    fname_str = sprintf('stability_sim_from_pred_trav_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d',...
           num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));       
    fname = fullfile(load_folder, fname_str);
    load(fname, 'trav_wave_all', 'trav_wave_all_strict', 'unmodified_all');
    
    % get list of simulations that generated trav. waves
    trav_wave_conds_met_sim = (trav_wave_all & unmodified_all);
    fprintf('Network %d, #trav. wave conditions = %d \n', network, sum(trav_wave_conds_met_sim) )
    
    Con_wave_sim = squeeze( Con_all_network(trav_wave_conds_met_sim, :) );
    K_wave_sim = squeeze(K_all_network(trav_wave_conds_met_sim, :, :));    
    
    % get K, Con parameters found by both theory and simulations
    idx_sel = squeeze(trav_wave_conds_met_theory & trav_wave_conds_met_sim);
    Con_wave_both = Con_all_network( idx_sel, :);
    K_wave_both = K_all_network(idx_sel, :, :);  
    
    % get K, Con parameters false positives (theory: TW, simulation: not
    % TW)
    idx_sel = squeeze(trav_wave_conds_met_theory & ~trav_wave_conds_met_sim);
    Con_wave_FP = Con_all_network( idx_sel, :);
    K_wave_FP = K_all_network(idx_sel, :, :);  
    
    % get K, Con parameters false negatives (theory: Not TW, simulation:
    % TW)
    idx_sel = squeeze(~trav_wave_conds_met_theory & trav_wave_conds_met_sim);
    Con_wave_FN = Con_all_network( idx_sel, :);
    K_wave_FN = K_all_network(idx_sel, :, :);  
%end
%}

%% ----------Study effect of varying parameters----------------------------
% Fix pw
pw = [2/gz 2/gz];

% Vary gz
gz_all = 5:5:50; %[5 15 50]; %[5 15 50]; %[5 15 50]; %5:5:50;
% Vary a0
a0_all = 0.1:0.1:4; %[0.5 1.5 5]; %0.1:0.01:5; %0.1:0.1:5;
% Vary lambda1, lambda2
%lambda1_all = 1:0.1:2; 
%lambda2_all = 1:0.1:5; 

l1_all = 0.1:0.1:6;
l2_all = 0.1:0.1:4; %4; %:5;

% results to store
Q_value_all = zeros(numel(gz_all), numel(a0_all), numel(l1_all), numel(l2_all));
Q_value_N_inf_all = zeros(numel(gz_all), numel(a0_all), numel(l1_all), numel(l2_all));
Q_value_L_inf_all = zeros(numel(gz_all), numel(a0_all), numel(l1_all), numel(l2_all));
Q_value_N_inf_L_inf_all = zeros(numel(gz_all), numel(a0_all), numel(l1_all), numel(l2_all));
Q_value_NNA_all = zeros(numel(gz_all), numel(a0_all), numel(l1_all), numel(l2_all));

for j=1:numel(gz_all)
    this_gz = gz_all(j);
    disp(this_gz);
    this_pw = [1/8 1/8]; %[2/this_gz 2/this_gz];
    
    % Calculate distances
    [~, this_dist] = initial_cells_random_markov_periodic(this_gz, mcsteps, rcell);
    
    for i=1:numel(a0_all)
        this_a0 = a0_all(i);
        disp(this_a0);
        for i1=1:numel(l1_all)
            for i2=1:numel(l2_all)
                %this_lambda = [lambda1_all(i1) lambda2_all(i2)];
                this_l = [l1_all(i1) l2_all(i2)];
                this_lambda = this_l*this_a0;
                % -----------Calculate interaction strengths-----------------------
                [fN, fnn, ~] = calc_fN(this_a0, rcell, this_gz, this_dist, this_lambda);
                Rcell = rcell*this_a0;
                f_MF = (4*pi/this_a0^2/sqrt(3)).*exp(Rcell./this_lambda).*sinh(Rcell./this_lambda).*exp(-3*this_a0./(2*this_lambda));
                %------------------------------------------------------------------

                % (1) Fixed phase space region [1,L]x[1,L]
                L = 1000;
                A_frac_self = calc_area(L, fN, fnn, pw, 'self', 1, 0, f_MF);
                A_frac_mutual = zeros(2, 1);
                A_frac_mutual(1) = calc_area(L, fN, fnn, pw, 'mutual', 2, 0, f_MF);
                A_frac_mutual(2) = calc_area(L, fN, fnn, pw, 'mutual', 1, 0, f_MF);
                Q_value = A_frac_self*prod(A_frac_mutual);
                % disp(Q_value);
                %{
                % (2) Limit N -> infty
                L = 1000;
                A_frac_self = calc_area(L, fN, fnn, pw, 'self', 1, 1, f_MF);
                A_frac_mutual = zeros(2, 1);
                A_frac_mutual(1) = calc_area(L, fN, fnn, pw, 'mutual', 2, 1, f_MF);
                A_frac_mutual(2) = calc_area(L, fN, fnn, pw, 'mutual', 1, 1, f_MF);
                Q_value_N_inf = A_frac_self*prod(A_frac_mutual);
                % disp(Q_value);

                % (3) limit L -> infty
                A_frac_self = calc_area_inf(L, fN, fnn, pw, 'self', 1, 0, f_MF);
                A_frac_mutual = zeros(2, 1);
                A_frac_mutual(1) = calc_area_inf(L, fN, fnn, pw, 'mutual', 2, 0, f_MF);
                A_frac_mutual(2) = calc_area_inf(L, fN, fnn, pw, 'mutual', 1, 0, f_MF);
                Q_value_L_inf = A_frac_self*prod(A_frac_mutual);
                % disp(Q_value);

                % (4) limit L -> infty, N -> infty
                A_frac_self = calc_area_inf(L, fN, fnn, pw, 'self', 1, 1, f_MF);
                A_frac_mutual = zeros(2, 1);
                A_frac_mutual(1) = calc_area_inf(L, fN, fnn, pw, 'mutual', 2, 1, f_MF);
                A_frac_mutual(2) = calc_area_inf(L, fN, fnn, pw, 'mutual', 1, 1, f_MF);
                Q_value_L_inf_N_inf = A_frac_self*prod(A_frac_mutual);

                % (5) limit L -> infty, N -> infty, NNA
                f_MF = 0;
                A_frac_self = calc_area_inf(L, fN, fnn, pw, 'self', 1, 1, f_MF);
                A_frac_mutual = zeros(2, 1);
                A_frac_mutual(1) = calc_area_inf(L, fN, fnn, pw, 'mutual', 2, 1, f_MF);
                A_frac_mutual(2) = calc_area_inf(L, fN, fnn, pw, 'mutual', 1, 1, f_MF);
                Q_NNA = A_frac_self*prod(A_frac_mutual);
                %}
                %------------------------------
                % Store results
                Q_value_all(j, i, i1, i2) = Q_value;
                %{
                Q_value_N_inf_all(j, i, i1, i2) = Q_value_N_inf;
                Q_value_L_inf_all(j, i, i1, i2) = Q_value_L_inf;
                Q_value_N_inf_L_inf_all(j, i, i1, i2) = Q_value_L_inf_N_inf;
                Q_value_NNA_all(j, i, i1, i2) = Q_NNA;
                %}
            end
        end
    end
end

% Save data 
fname_str = strrep(sprintf('Q_values_gz_%d_to_%d_a0_%.1f_to_%.1f_l1_%.1f_to_%.1f_l2_%.1f_to_%.1f',...
    gz_all(1), gz_all(end), a0_all(1), a0_all(end), l1_all(1), l1_all(end), l2_all(1), l2_all(end)), '.', 'p');
fname = fullfile( save_folder, 'data', fname_str );
save( fname, 'gz_all', 'a0_all', 'l1_all', 'l2_all', 'Q_value_all' );

%% Plot and save examples of phase plots
% Fix pw
pw = [1/8 1/8];
L = 1000;
this_gz = 15;
this_a0 = 0.2;
this_lambda = lambda;

[~, this_dist] = initial_cells_random_markov_periodic(this_gz, mcsteps, rcell);
[fN, fnn, ~] = calc_fN(this_a0, rcell, this_gz, this_dist, this_lambda);
showlines = 0;

h = plot_boundaries_K_Con(fN, fnn, pw, L, showlines);

% Save figure
qsave = 1;
fname_str = strrep(sprintf('Phase_boundaries_gz_%d_a0_%.1f_lambda_%.1f_%.1f_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d',...
	this_gz, this_a0, this_lambda(1), this_lambda(2), ...
    num_waves, wave_type, network, states_perm(1), states_perm(2),...
	states_perm(3), states_perm(4)), '.', 'p');
fname = fullfile(save_folder, 'phase_boundaries', fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%% Plot Q-value vs a0
idx_gz = 1; % fix gz
this_gz = gz_all(idx_gz);
idx_l = [1 1]; % fix lambda
this_l = [l1_all(idx_l(1)) l2_all(idx_l(2))];

h=figure;
hold on
ls = '--';
clrs = {'b', 'r', 'g', 'm'};
p1=plot(a0_all, Q_value_all(idx_gz, :, idx_l(1), idx_l(2)), '-', 'LineWidth', 3);
p2=plot(a0_all, Q_value_N_inf_all(idx_gz, :, idx_l(1), idx_l(2)), ls, 'LineWidth', 2);
p3=plot(a0_all, Q_value_L_inf_all(idx_gz, :, idx_l(1), idx_l(2)), ls, 'LineWidth', 2);
p4=plot(a0_all, Q_value_N_inf_L_inf_all(idx_gz, :, idx_l(1), idx_l(2)), ls, 'LineWidth', 2);
p5=plot(a0_all, Q_value_NNA_all(idx_gz, :, idx_l(1), idx_l(2)), ls, 'LineWidth', 2);
xlabel('a_0');
ylabel('Q-value');
legend([p1 p2 p3 p4 p5], {sprintf('L=%d, N=%d', L, this_gz^2), sprintf('L=%d, N=\\infty', L),...
     sprintf('L=\\infty, N=%d', this_gz^2) 'L=\infty, N=\infty', 'L=\infty, N=\infty, NNA'});
set(gca, 'FontSize', 32);
xlim([0 max(a0_all)]);

% Save figure
qsave = 0;
fname_str = sprintf('Q_value_vs_a0_gz_%d_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d',...
	this_gz, num_waves, wave_type, network, states_perm(1), states_perm(2),...
	states_perm(3), states_perm(4));
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%% Plot Q-value vs gz
idx_a0 = 1; % fix a0
this_a0 = a0_all(idx_a0);
idx_l = [1 1]; % fix lambda
this_l = [l1_all(idx_l(1)) l2_all(idx_l(2))];

h=figure;
hold on
plot(gz_all, Q_value_all(:, idx_a0, idx_l(1), idx_l(2)), 'o--', 'LineWidth', 1.5)
plot(gz_all, Q_value_L_inf_all(:, idx_a0, idx_l(1), idx_l(2)), 'o--', 'LineWidth', 1.5)
plot(gz_all, Q_value_N_inf_L_inf_all(:, idx_a0, idx_l(1), idx_l(2)), '--', 'LineWidth', 1.5)
plot(gz_all, Q_value_NNA_all(:, idx_a0, idx_l(1), idx_l(2)), '--', 'LineWidth', 1.5)
xlabel('Grid size');
ylabel('Q-value');
set(gca, 'FontSize', 32);
legend({sprintf('L=%d', L), 'L=\infty',...
    'L=\infty, N=\infty', 'L=\infty, N=\infty, NNA'}, 'Location', 'ne');
ylim([0 0.05]);

% Save figure
qsave = 0;
fname_str = strrep(sprintf('Q_value_vs_gz_a0_%.1f_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d',...
	this_a0, num_waves, wave_type, network, states_perm(1), states_perm(2),...
	states_perm(3), states_perm(4)), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot (fractional) errors
% -> correct indices
idx_gz = 3;
this_gz = gz_all(idx_gz);

frac_err_N_inf = abs(Q_value_N_inf_all-Q_value_all)./Q_value_all;
frac_err_L_inf = abs(Q_value_L_inf_all-Q_value_all)./Q_value_all;
frac_err_L_inf_N_inf = abs(Q_value_N_inf_L_inf_all-Q_value_all)./Q_value_all;
frac_err_NNA = abs(Q_value_NNA_all-Q_value_all)./Q_value_all;
ls = '--';

h=figure;
hold on
plot([0 0], [0 0]);
p2 = plot(a0_all(1:end), frac_err_N_inf(idx_gz, 1:end), ls, 'LineWidth', 2 );
p3 = plot(a0_all(1:end), frac_err_L_inf(idx_gz, 1:end), ls, 'LineWidth', 2 );
p4 = plot(a0_all(1:end), frac_err_L_inf_N_inf(idx_gz, 1:end), ls, 'LineWidth', 2 );
p5 = plot(a0_all(1:end), frac_err_NNA(idx_gz, 1:end), ls, 'LineWidth', 2 );
legend([p2 p3 p4 p5], {sprintf('L=%d, N=\\infty', L), sprintf('L=\\infty, N=%d', this_gz^2),...
    'L=\infty, N=\infty', 'L=\infty, N=\infty, NNA'});
xlabel('a_0');
ylabel('Fractional error');
set(gca, 'FontSize', 32);
ylim([0 1]);

% Save figure
qsave = 0;
fname_str = strrep(sprintf('Q_value_frac_error_vs_a0_gz_%d_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d',...
	this_gz, num_waves, wave_type, network, states_perm(1), states_perm(2),...
	states_perm(3), states_perm(4)), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);


%% Plot Q-value vs gz, a0 simultaneously 


%% Plot Q-value vs lambda
idx_a0 = 1; % fix a0
idx_gz = 1; % fix gz
this_a0 = a0_all(idx_a0);
this_gz = gz_all(idx_gz);

h=figure;
hold on
Q_data = squeeze(Q_value_all(idx_gz, idx_a0, :, 1:end));
imagesc(l2_all, l1_all(1:end), Q_data)
set(gca, 'YDir', 'normal')
%plot(gz_all, Q_value_L_inf_all(idx_gz, idx_a0, :, :), 'o--', 'LineWidth', 1.5)
%plot(gz_all, Q_value_N_inf_L_inf_all(idx_gz, idx_a0, :, :), '--', 'LineWidth', 1.5)
%plot(gz_all, Q_value_NNA_all(idx_gz, idx_a0, :, :), '--', 'LineWidth', 1.5)
xlabel('l^{(2)} = \lambda^{(2)}/a_0');
ylabel('l^{(1)} = \lambda^{(1)}/a_0');
set(gca, 'FontSize', 32);
%title(sprintf('gz = %d, a0 = %.2f', this_gz, this_a0));
title(sprintf('gz = %d', this_gz));
%legend({sprintf('L=%d', L), 'L=\infty',...
%    'L=\infty, N=\infty', 'L=\infty, N=\infty, NNA'}, 'Location', 'ne');
%ylim([0 0.05]);
colorbar;

% Find maximum
[i1, i2] = find(Q_data == max(Q_data(:)));
fprintf('max(Q) = %.3f, l1 = %.2f, l2 = %.2f \n', max(Q_data(:)), l1_all(i1), l2_all(i2) )
% Plot maximum
scatter(  l2_all(i2), l1_all(i1), 'rx' );

% Save figure
qsave = 0;
fname_str = strrep(sprintf(...
    'Q_value_vs_lambda_gz_%d_network_%d_range1_%.2f_%.2f_range1_%.2f_%.2f',...
	this_gz, network, l1_all(1), l1_all(end), l2_all(1), l2_all(end)), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Optimal robustness vs a0
% Maximize robustness function L->infty, N->infty
gz_all = 5:5:50;
a0_all = 0.2:0.01:1.5;
a0_max_all = zeros(numel(gz_all), 1);
a0_max_NNA_all = zeros(numel(gz_all), 1);
Q_value_max_all = zeros(numel(gz_all), 1);
Q_value_max_NNA_all = zeros(numel(gz_all), 1);

for j=1:numel(gz_all)
    gz = gz_all(j);
    disp(gz);
    %this_pw = [2/this_gz 2/this_gz];
    
    % Calculate distances
    [~, this_dist] = initial_cells_random_markov_periodic(this_gz, mcsteps, rcell);

    Q_value_all_temp = zeros(numel(a0_all), 1);
    Q_value_NNA_all_temp = zeros(numel(a0_all), 1);
    for i=1:numel(a0_all)
        this_a0 = a0_all(i);

        % -----------Calculate interaction strengths-----------------------
        [fN, fnn, ~] = calc_fN(this_a0, rcell, this_gz, this_dist, lambda);
        Rcell = rcell*this_a0;
        f_MF = (4*pi/this_a0^2/sqrt(3)).*exp(Rcell./lambda).*sinh(Rcell./lambda).*exp(-3*this_a0./(2*lambda));
        % -> check dependence on lambda
        %------------------------------------------------------------------
        % (1) finite region, finite N
        L = 1000;
        A_frac_self = calc_area(L, fN, fnn, pw, 'self', 1, 0, f_MF);
        A_frac_mutual = zeros(2, 1);
        A_frac_mutual(1) = calc_area(L, fN, fnn, pw, 'mutual', 2, 0, f_MF);
        A_frac_mutual(2) = calc_area(L, fN, fnn, pw, 'mutual', 1, 0, f_MF);
        Q_value = A_frac_self*prod(A_frac_mutual);
        % disp(Q_value);

        % (5) limit L -> infty, N -> infty, NNA
        %{
        f_MF = 0;
        A_frac_self = calc_area_inf(L, fN, fnn, pw, 'self', 1, 1, f_MF);
        A_frac_mutual = zeros(2, 1);
        A_frac_mutual(1) = calc_area_inf(L, fN, fnn, pw, 'mutual', 2, 1, f_MF);
        A_frac_mutual(2) = calc_area_inf(L, fN, fnn, pw, 'mutual', 1, 1, f_MF);
        Q_NNA = A_frac_self*prod(A_frac_mutual);
        %}

        Q_value_all_temp(i) = Q_value;
        %Q_value_NNA_all_temp(i) = Q_NNA;
    end
    
    % Store results
    a0_max_all(j) = a0_all( Q_value_all_temp == max(Q_value_all_temp) );
    Q_value_max_all(j) = Q_value_all_temp( Q_value_all_temp == max(Q_value_all_temp) );

    %a0_max_NNA_all(j) = a0_all( Q_value_NNA_all_temp == max(Q_value_NNA_all_temp) );
    %Q_value_max_NNA_all(j) = Q_value_NNA_all_temp( Q_value_NNA_all_temp == max(Q_value_NNA_all_temp) );
end

%% Plot solution
%{
figure;
plot(a0_all, Q_value_all);
fprintf('argmax(Q-value) = %.3f \n', a0_max);
%}

% Plot a0 at max Q
h1=figure;
box on
hold on
plot(gz_all, a0_max_all, 'o-', 'LineWidth', 2);
%plot(gz_all, a0_max_NNA_all);
xlabel('Grid Size');
ylabel('Optimal a_0');
set(gca, 'FontSize', 32);
ylim([0 5]);

% Plot max Q
h2=figure;
box on
plot(gz_all, Q_value_max_all, 'o-', 'LineWidth', 2);
%plot(gz_all, Q_value_NNA_max_all);
xlabel('Grid Size');
ylabel('Max Q-value');
set(gca, 'FontSize', 32);
ylim([0 0.05]);

%% print results
%fprintf('%d argmax(Q-value) = %.3f \n', a0_max);
t = table(gz_all', Q_value_max_all, a0_max_all, 'VariableNames', {'gz', 'Qvalue', 'a0'});

%%

% Save figures
qsave = 1;
fname_str = strrep(sprintf('Optimal_Q_value_argmax_a0_vs_gz_%d_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d',...
	this_gz, num_waves, wave_type, network, states_perm(1), states_perm(2),...
	states_perm(3), states_perm(4)), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h1, 10, 8, fname, '.pdf', qsave);

qsave = 1;
fname_str = strrep(sprintf('Optimal_Q_value_vs_gz_%d_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d',...
	this_gz, num_waves, wave_type, network, states_perm(1), states_perm(2),...
	states_perm(3), states_perm(4)), '.', 'p');
fname = fullfile(save_folder, fname_str);
save_figure(h2, 10, 8, fname, '.pdf', qsave);
%% Troubleshooting
this_l = [2 1];
this_lambda = this_l*this_a0;
% -----------Calculate interaction strengths-----------------------
[fN, fnn, ~] = calc_fN(this_a0, rcell, this_gz, this_dist, this_lambda);
Rcell = rcell*this_a0;
f_MF = (4*pi/this_a0^2/sqrt(3)).*exp(Rcell./this_lambda).*sinh(Rcell./this_lambda).*exp(-3*this_a0./(2*this_lambda));
%------------------------------------------------------------------

[A, C1, C2] = calc_area_inf(L, fN, fnn, pw, interaction, molecule, N_inf, f_MF);

% (1) Fixed phase space region [1,L]x[1,L]
%{
L = 1000;
A_frac_self = calc_area(L, fN, fnn, pw, 'self', 1, 0, f_MF);
A_frac_mutual = zeros(2, 1);
A_frac_mutual(1) = calc_area(L, fN, fnn, pw, 'mutual', 2, 0, f_MF);
A_frac_mutual(2) = calc_area(L, fN, fnn, pw, 'mutual', 1, 0, f_MF);
Q_value = A_frac_self*prod(A_frac_mutual);

disp(Q_value)
%}
%{
L = 1000;
a0 = 0.15;
[fN, fnn, ~] = calc_fN(a0, rcell, gz, dist, lambda);

%disp( calc_area(L, fN, fnn, pw, 1, 1)/((L-1)^2) )
%disp( calc_area(L, fN, fnn, pw, 1, 2)/((L-1)^2) )
%disp( calc_area(L, fN, fnn, pw, 2, 1)/((L-1)^2) )
% 1<-2 and 2<-1 interactions
C1 = (L - 2*fnn - (fN'-6*fnn).*(1-2*pw))./(1+4*fnn + (fN'-6*fnn).*pw);
I1 = (C1.^2/2 - C1 + 1/2).*(1+2*fnn); % integral 1
I2 = L*(L-C1) - ( (2*fnn + (fN'-6*fnn).*pw)/2.*(L.^2-C1.^2) +...
    (1+4*fnn+(fN'-6*fnn).*(1-pw)).*(L-C1) ) ; % integral 2 -> check calculation
V_frac_mutual = (I1 + I2)/(L-1)^2; 
% 1<-1 interaction
V_frac_self = (fnn(1).*(L^2 - 2*L) - fnn(1))./((L-1)^2);

Rcell = rcell*a0;
f_MF = (4*pi/a0^2/sqrt(3)).*exp(Rcell./lambda).*sinh(Rcell).*exp(-3*a0./(2*lambda));
disp([(fN'-6*fnn); f_MF]);

disp([...%V_frac_self ...
    calc_area(L, fN, fnn, pw, 'self', 1, 0, f_MF) ...
    calc_area(L, fN, fnn, pw, 'self', 1, 1, f_MF) ...
    calc_area_inf(L, fN, fnn, pw, 'self', 1, 0, f_MF) ...
    calc_area_inf(L, fN, fnn, pw, 'self', 1, 1, f_MF) ...
    calc_area_inf(L, fN, fnn, pw, 'self', 1, 1, 0)]);
disp([...%V_frac_mutual(2) ...
    calc_area(L, fN, fnn, pw, 'mutual', 2, 0, f_MF) ...
    calc_area(L, fN, fnn, pw, 'mutual', 2, 1, f_MF) ...
    calc_area_inf(L, fN, fnn, pw, 'mutual', 2, 0, f_MF) ...
    calc_area_inf(L, fN, fnn, pw, 'mutual', 2, 1, f_MF) ...
    calc_area_inf(L, fN, fnn, pw, 'mutual', 2, 1, 0)]);
disp([...%V_frac_mutual(1) ...
    calc_area(L, fN, fnn, pw, 'mutual', 1, 0, f_MF) ...
    calc_area(L, fN, fnn, pw, 'mutual', 1, 1, f_MF) ...
    calc_area_inf(L, fN, fnn, pw, 'mutual', 1 ,0, f_MF) ...
    calc_area_inf(L, fN, fnn, pw, 'mutual', 1 ,1, f_MF) ...
    calc_area_inf(L, fN, fnn, pw, 'mutual', 1 ,1, 0)]);
        
%disp(Q_value)
disp('--')
%}
%% Functions
function [A, C1, C2] = calc_area_inf(L, fN, fnn, pw, interaction, molecule, N_inf, f_MF)
    % Calculate in the L->infty limit
    % N_inf = 1: N->infty limit
    % interaction: 'self' -> self (1); 'mutual' -> mutual (2)
    % molecule: 1 or 2
    
    % convert string input to numeric
    int = find( cellfun(@(x) ~isempty(x), regexp( interaction, {'self', 'mutual'}, 'match' )), 1);
    
    % specify molecule 
    fN = fN(molecule);
    fnn = fnn(molecule);
    pw = pw(molecule);
    
    if ~N_inf % if not N-> infty
        f_MF = (fN'-6.*fnn); % restore original expression
    end
    
    % calculate intersection points and define integrals
    % integrals are normalized by (L-1)^2, arguments are alpha factors (for L, alpha=1)
    syms x x1 x2 % for integrals
    switch int
        case 1 % self
            % intersection points 
            C1 = (L - 1 - 4.*fnn - (fN' - 6*fnn).*(1-pw))./(2.*fnn+f_MF.*pw);
            C2 = (L - 1 - 6.*fnn - (fN' - 6*fnn).*(1-pw))./(f_MF.*pw);
            
            % alpha values
            alpha1 = 1./(2.*fnn+f_MF.*pw);
            alpha2 = 1./(f_MF.*pw);
            
            % Integrals
            I1 = @(a) 1/2.*(1+2*fnn)*a^2; % integral 1
            I2 = @(a1, a2) (a2-a1) - (fnn + f_MF.*pw/2).*(a2.^2-a1.^2);
        case 2 % mutual
            % intersection points 
            C1 = (L - 0 - 2.*fnn - (fN' - 6*fnn).*(1-pw))./(1 + 4.*fnn+f_MF.*pw);
            C2 = (L - 1 - 6.*fnn - (fN' - 6*fnn).*(1-pw))./(2*fnn + f_MF.*pw);
            
            % alpha values
            alpha1 = 1./(1 + 4.*fnn+f_MF.*pw);
            alpha2 = 1./(2*fnn + f_MF.*pw);
            
            % Integrals
            I1 = @(a) a^2.*fnn;
            I2 = @(a1, a2) a2-a1 - f_MF.*pw/2.*(a2.^2-a1.^2);
    end
    
    % Cases
    if C1<L && C2<L
        % case A
        %disp('case A');
        A = I1(alpha1) + I2(alpha1, alpha2);
    elseif C1<L && C2>L
        % case B   
        %disp('case B');
        A = I1(alpha1) + I2(alpha1, 1);
    else
        % case C
        %disp('case C');
        A = I1(1);
    end
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
    %{
    % intersection points 
    C1_self = (L - 1 - 4.*fnn - f_MF.*(1-pw))./(2.*fnn+f_MF.*pw);
    C2_self = (L - 1 - 6.*fnn - f_MF.*(1-pw))./(f_MF.*pw);
    C1_mutual = (L - 0 - 2.*fnn - f_MF.*(1-pw))./(1 + 4.*fnn+f_MF.*pw);
    C2_mutual = (L - 1 - 6.*fnn - f_MF.*(1-pw))./(2*fnn + f_MF.*pw);
    C1_both = {C1_self, C1_mutual};
    C2_both = {C2_self, C2_mutual};
    C1 = C1_both{int}(molecule);
    C2 = C2_both{int}(molecule);
    
    % Integrals
    syms x x1 x2
    I1_mutual = @(x) (x.^2 - 2*x + 1)/2.*(1+2*fnn); % integral 1
    I2_mutual = @(x1, x2) L*(x2-x1) - ( (fnn + f_MF.*pw/2).*(x2.^2-x1.^2) -...
        (1 + 4*fnn + f_MF.*(1-pw)).*(x2-x1) ) ; % integral 2 -> check calculation
    I1_self = @(x) (x.^2 - 2*x - 1).*fnn;
    I2_self = @(x1, x2) L*(x2-x1) - f_MF.*pw/2.*(x2.^2-x1.^2) -...
        (1 + 6*fnn + f_MF.*(1-pw)).*(x2-x1);
    I1 = {I1_self, I1_mutual};
    I2 = {I2_self, I2_mutual};

    % Cases
    if C1<L && C2<L
        % case A
        %disp('case A');
        A = I1{int}(C1) + I2{int}(C1, C2);
    elseif C1<L && C2>L
        % case B   
        %disp('case B');
        A = I1{int}(C1) + I2{int}(C1, L);
    else
        % case C
        %disp('case C');
        A = I1{int}(L);
    end
    A = A(molecule)/(L-1)^2;
    %}
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
end

function [M_int_found] = get_found_M_int(y_found)
    % Get a list of all interaction matrices for the found networks
    M = [0 1 -1]; % index to interaction
    M_int_all_reduced = {};
    done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
    for k=1:3^4
        [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], k);
        gM = [i22 i21; i12 i11];
        M_int = [M(i11) M(i12); M(i21) M(i22)];
        if done(i11,i12,i21,i22)
            continue
        elseif k==1
            continue
        else
            M_int_all_reduced{end+1} = M_int;
            done(i11,i12,i21,i22) = 1;
            done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
        end
    end

    M_int_found = cell(numel(y_found), 1);
    for i=1:numel(y_found)
        M_int_found{i} = M_int_all_reduced{y_found(i)};
        %disp(M_int_found{i})
    end
end

function h = plot_Con_K_map_minimal(idx_gene, idx_mol)
    % Plot phase diagram without colors
    % idx_gene: index of gene under control
    % idx_mol: index of sensed molecule / molecule affecting the gene (1 <= i <= L)
        
    % Parameter range map
    Con_max = 1000;
    K_max = 1000;

    % plot figure
    h = figure;
    hold on
    %set(him, 'AlphaData', out > 0); % invisible if not from any region
    
    % adjust the graph
    set(gca,'ydir', 'normal', 'FontSize', 24)
    xlabel(sprintf('K^{(%d%d)}', idx_gene, idx_mol), 'FontSize', 24)
    ylabel(sprintf('C_{ON}^{(%d)}', idx_mol), 'FontSize', 24)
    ylim([1 Con_max])
    xlim([1 K_max])
    %title(sprintf('$$f_N = %.3f, a_0 = %.2f$$', fN, a0), 'FontSize', 30)
    title(sprintf('Interaction %d \\leftarrow %d', idx_gene, idx_mol))
    %title(sprintf('Molecule %d', idx_mol));    
end

function h = plot_phase_diagram_local(a0, rcell, lambda, dist, idx_gene, idx_mol)
    % idx_gene: index of gene under control
    % idx_mol: index of sensed molecule / molecule affecting the gene (1 <= i <= L)
    
    % calculate fN
    Rcell = rcell*a0;
    fN = zeros(2,1);
    %[dist, ~] = init_dist_hex(gz, gz);
    dist_vec = a0*dist(1,:);
    r = dist_vec(dist_vec>0); % exclude self influence
    for i=1:2 % calculate signaling strength
        fN(i) = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(i)).*(lambda(i)./r)) ); 
    end

    % Parameter range map
    num_points = 1000;
    Con_max = 1000;
    K_max = 1000;
    Con_vec = linspace(1, Con_max, num_points);
    K_vec = linspace(1, K_max, num_points);
    [K, Con] = meshgrid(K_vec, Con_vec);

    % Make 4 limiting regions as boolean matrices
    R1 = (1+fN(idx_mol) - K) > 0; % Everything ON
    R2 = ((1+fN(idx_mol))*Con - K ) < 0; % Everything OFF
    R3 = ((Con + fN(idx_mol) - K) > 0 & (1+fN(idx_mol) - K) < 0); % ON remains ON & not all ON
    R4 = ((1+ fN(idx_mol)*Con - K) < 0 & ((1+fN(idx_mol))*Con - K ) > 0) ; % OFF remains OFF & not all OFF
    %R3 = (Con > K & (K-1)./Con > fN); % autonomous cells for Son > K
    %R4 = (Con <= K & K - Con < fN & (K-1)./Con > fN); % autonomous cells for Son < K

    out = R1 + 2*R2 + 3*R3 + 4*R4; % only regions 3 and 4 can overlap
    if ~isempty(find(unique(out)==0, 1))
        map_idx = 5; % activation-deactivation
        out(out==0) = 5; 
        %phase = 'non-A';
        phase = 'U';
    elseif ~isempty(find(unique(out)==7, 1)) 
        map_idx = 6; % autonomy
        out(out==7) = 5; % ON remains ON & OFF remains OFF
        phase = 'A01';
    end

    % plot figure
    h = figure;
    hold on
    him = imagesc(K_vec, Con_vec, out);

    %set(him, 'AlphaData', out > 0); % invisible if not from any region
    % R1 -> black
    % R2 -> white
    % R3 -> green
    % R4 -> red
    % activation-deactivation -> magenta
    % autonomy -> gray
    map = [0, 0, 0
        1, 1, 1
        0, 1, 0
        1, 0, 0
        1, 1, 0
        0.5, 0.5, 0.5];
    tmp = map([1:4 map_idx], :);
    colormap(tmp);
    c=colorbar;
    set(c, 'YTick', 1+2/5+4/5*(0:4));

    %phase_labels = {'all>K','all<K','A1','A0',phase};
    phase_labels = {'P1','P0','A1','A0',phase};
    set(c, 'TickLabels', phase_labels);

    % adjust the graph
    set(gca,'ydir', 'normal', 'FontSize', 24)
    xlabel('K', 'FontSize', 24)
    ylabel('$$C_{ON}$$', 'FontSize', 24)
    ylim([1 Con_max])
    xlim([1 K_max])
    %title(sprintf('$$f_N = %.3f, a_0 = %.2f$$', fN, a0), 'FontSize', 30)
    title(sprintf('Interaction $$%d \\leftarrow %d$$', idx_gene, idx_mol))
    %title(sprintf('Molecule %d', idx_mol));
end

%% Previous code

% -------------------------------------------------------------------------
% Calculation that does not properly take into account the various
% possibilities for intersecting the boundaries 
% -------------------------------------------------------------------------

%{
%% (1) Calculate phase space region for given region of size L
L = 1000;

% 1<-2 and 2<-1 interactions
C1 = (L - 2*fnn - (fN'-6*fnn).*(1-2*pw))./(1+4*fnn + (fN'-6*fnn).*pw);
I1 = (C1.^2/2 - C1 + 1/2).*(1+2*fnn); % integral 1
I2 = L*(L-C1) - ( (2*fnn + (fN'-6*fnn).*pw)/2.*(L.^2-C1.^2) +...
    (1+4*fnn+(fN'-6*fnn).*(1-pw)).*(L-C1) ) ; % integral 2 -> check calculation
V_frac_mutual = (I1 + I2)/(L-1)^2; 
% 1<-1 interaction
V_frac_self = (fnn(1).*(L^2 - 2*L) - fnn(1))./((L-1)^2);
%disp(V_frac_self);
%disp(V_frac_mutual(2));
%disp(V_frac_mutual(1));

Q_value = prod(V_frac_mutual)*V_frac_self;
disp(Q_value);
%% (2) Calculate phase space region in limit L -> infty
% (1) rectangle with width Km = upper bound K at given Cm
%{
V12 = (1/2+fnn(2))/(1+4*fnn(2)+(fN(2) - 6*fnn(2))*pw(2));
V21 = V12;
V22 = fnn(2)/(2*fnn(2) + (fN(2) - 6*fnn(2))*pw(2));
disp(V12*V21*V22)
%}

% (2) rectangle [1, L] x [1, L]
alpha = 1./(1+4*fnn + (fN'-6*fnn).*pw);
V_frac_mutual = (1/2+fnn).*alpha.^2 + (1-alpha) - (fnn+(fN'-6*fnn).*pw/2).*(1-alpha.^2);
V_frac_self = fnn(1); %/(2*fnn(2) + (fN(2) - 6*fnn(2))*pw(2));
Q_value_L_inf = prod(V_frac_mutual)*V_frac_self;
disp(Q_value_L_inf)

%% (3) Calculate phase space region in limit L -> infty, N -> infty
Rcell = rcell*a0;
f_MF = (4*pi/a0^2/sqrt(3)).*exp(Rcell./lambda).*sinh(Rcell).*exp(-3*a0./(2*lambda));

% Check calculation
%{
disp(f_MF)
A_cell = sqrt(3)/2*a0^2;
G = 2*pi*exp(Rcell./lambda)*sinh(Rcell);
disp( G/A_cell.*exp(-3*a0./(2*lambda)) )
disp( (fN'-6*fnn) )
L = a0*sqrt((N-7)*sqrt(3)/(2*pi) + 9/4); disp(exp(-L./lambda)); % -> negligibly small
%}

alpha = 1./(1+4*fnn + f_MF.*pw);
V_frac_mutual = (1/2+fnn).*alpha.^2 + (1-alpha) - (fnn+f_MF.*pw/2).*(1-alpha.^2);
V_frac_self = fnn(1); %/(2*fnn(2) + (fN(2) - 6*fnn(2))*pw(2));
Q_value_L_inf_N_inf = prod(V_frac_mutual)*V_frac_self;
disp(Q_value_L_inf_N_inf)

%% (4) Calculate phase space region in the nearest neighbour approximation (NNA)
% limit L -> infty, N -> infty
alpha_NNA = 1./(1+4*fnn);
Q_NNA = prod(-alpha_NNA/2 + 1 - fnn)*fnn(1);
disp(Q_NNA);
%}

%{
% Calculate fN vs a0
this_a0 = 1.5;
gz_all = 5:5:50;
fN_all = zeros(numel(gz_all), 2);
for i=1:numel(gz_all)
    this_gz = gz_all(i);
    [~, this_dist] = initial_cells_random_markov_periodic(this_gz, mcsteps, rcell);
    [fN, fnn, ~] = calc_fN(this_a0, rcell, this_gz, this_dist, lambda);
    fN_all(i, :) = fN;
end
h=figure;
plot(gz_all, fN_all);
%}