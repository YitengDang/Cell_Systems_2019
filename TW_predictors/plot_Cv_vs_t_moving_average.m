%% Plots the coefficient of variation (CV) over time for a moving average
clear all
close all
clc
set(0,'defaulttextinterpreter', 'tex');

%% Import data
%folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\filmstrip_selection';
%fname_str = 'travelling_wave_osc_background_network_20_broad_TW_vertical';
folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\Fig5-self-organisation';
load_fname_str = 'sample_complex_trajectory';

load(fullfile(folder, load_fname_str), 'cells_hist', 'positions', 'distances', 'save_consts_struct');
rcell = save_consts_struct.rcell;
a0 = save_consts_struct.a0;
t_out = numel(cells_hist)-1;

[period_ub, t_onset] = periodicity_test_short(cells_hist);
[period, t_onset] = periodicity_test_detailed(cells_hist, t_onset,...
    period_ub);
%% Get and plot p, I vs t
t0 = 0;
fig_pos = [1 1 10 8];

% plot p(t)
[msg, p_t] = plot_p_vs_t(cells_hist, t0, fig_pos);
disp(msg);

% plot I(t)
option = 1;
[msg, I_t] = plot_I_vs_t(cells_hist, t0, a0, distances, option, fig_pos);
disp(msg);

%% Calculate moving average
t_window = 10*period; % size of the window function
p_av = zeros(t_out-t_window, 2);
I_av = zeros(t_out-t_window, 2);
cv_p_av = zeros(t_out-t_window, 2); % CV for p(t)
cv_I_av = zeros(t_out-t_window, 2); % CV for I(t)

for tt=0:t_out-t_window+1
    p_av(tt+1, :) = mean(p_t(tt+1:tt+t_window, :), 1);
    I_av(tt+1, :) = mean(I_t(tt+1:tt+t_window, :), 1);

    cv_p_av(tt+1, :) = std(p_t(tt+1:tt+t_window, :), [], 1)./p_av(tt+1, :);
    cv_I_av(tt+1, :) = std(I_t(tt+1:tt+t_window, :), [], 1)./I_av(tt+1, :);
end

%% Plots (1) 
%{
% Plot averaged p_av(t), I_av(t)
h = figure;
plot(0:t_out-t_window+1, p_av)
ylim([0 1]);

h2 = figure;
plot(0:t_out-t_window+1, I_av)
ylim([-0.2 1]);
%}
%% Plots (2) 
save_folder = 'H:\My Documents\Multicellular automaton\figures\TW_predictors';

% Plot CV vs t
h3 = figure;
plot(0:t_out-t_window+1, cv_p_av, 'LineWidth', 2)
%ylim([0 1]);
xlabel('t');
ylabel('\sigma(p_t)/\mu(p_t)');
title(sprintf('t_{window} = %d', t_window));
set(gca, 'FontSize', 20);

fname_str = sprintf('cv_p_av_vs_t_%s_t_window_%d', load_fname_str, t_window);
fname = fullfile(save_folder, fname_str);
qsave = 1;
save_figure(h3, 10, 8, fname, '.pdf', qsave);
%%
h4 = figure;
plot(0:t_out-t_window+1, cv_I_av, 'LineWidth', 2)
%ylim([0 1]);
xlabel('t');
ylabel('\sigma(I_t)/\mu(I_t)');
title(sprintf('t_{window} = %d', t_window));
set(gca, 'FontSize', 20);

fname_str = sprintf('cv_I_av_vs_t_%s_t_window_%d', load_fname_str, t_window);
fname = fullfile(save_folder, fname_str);
qsave = 1;
save_figure(h4, 10, 8, fname, '.pdf', qsave);
