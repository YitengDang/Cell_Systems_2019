%% Plots the travelling wave abundances calculated from analytical estimate (Mathematica data)
clear all
close all
%% Load data
folder = 'H:\My Documents\Multicellular automaton\data\two_signals';
fname_str = 'NwaveDensity'; %'TravWaveDensity';
fname = fullfile(folder, fname_str);
xls_data = xlsread(fname);
%% Plot theoretical estimates vs N
idx_sel = 1:size(xls_data, 2);

h = figure;
plot(xls_data(1,idx_sel), 1./xls_data(2,idx_sel), 'bo-');
set(gca, 'YScale', 'log');
xlabel('Grid size');
%ylabel('Abundance');
ylabel('Time');
set(gca, 'FontSize', 24);

%% Parameters
idx_sel = 5;
gz = xls_data(1, idx_sel); % grid size

%% Plot cumulative probability for a Bernoulli process
% (Cumulative Poisson distribution)
digits(1000); % set precision
p = xls_data(2, idx_sel);

t_all = 10.^(48:0.1:54); %10.^(0:100);
prob_vs_t = (1-vpa(p)).^(t_all)*p;
q = 1-vpa(p);
prob_cumul_vs_t = (1-vpa(q).^t_all);

h=figure;
hold on
%plot(log10(t_all), prob_vs_t);
%plot(log10(t_all), prob_cumul_vs_t, 'o-');
plot(t_all, prob_cumul_vs_t, '-', 'LineWidth', 2);
plot([1./p 1./p], [0 1], '--');
%set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
ylim([0 1]);
%text(10^50, 0.7, sprintf('<t_{onset}>=%d', 1./p));
%xlabel('Time');
%ylabel('Cumulative probability');
set(gca, 'FontSize', 32);
box on

qsave = 1;
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_vs_N';
fname_str = sprintf('t_onset_TW_cumulative_prob_random_process_theory_gz%d_no_labels', gz);
fname = fullfile(folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%}
%% Plot simulation results
%gz_all = [8:2:20 25]; % from simulations

% New version that loads simulation data can be found in:
% -> analyze_trajectories_trav_wave_vs_variable 

%% Load processed data
%
load_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_vs_N';
fname_str = 'trav_wave_occur_vs_N_K12_10_nruns300_round_p_I_digits_5';
load(fullfile(load_path, strcat(fname_str, '.mat')), 't_onset_all');
%}
%% Plot distribution of t onset
t_data = t_onset_all(idx_sel, :);
edges = 0:max(t_data)+1;
hcounts = histcounts(t_data, edges, 'Normalization', 'cdf');
tmax = 6000;

h = figure;
hold on
plot(0:tmax, [hcounts ones(1, tmax-numel(hcounts)+1)] , 'LineWidth', 2);
plot([mean(t_data) mean(t_data)], [0 1], '--');
%xlabel('Time');
%ylabel('Cumulative probability');
set(gca, 'FontSize', 32);
set(gca, 'YTick', 0:0.2:1);
box on
% add mean as text
%text( 2000, 0.5, sprintf('<t_{onset}>=%d', mean(t_data) ) );

% save
qsave = 1;
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_vs_N';
fname_str = sprintf('t_onset_TW_cumulative_prob_simulations_gz%d_no_labels', gz);
fname = fullfile(folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);
