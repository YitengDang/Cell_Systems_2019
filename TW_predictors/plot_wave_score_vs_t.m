%% Plots one of the functions ("wave score") devised to study "waveness" of a simulation over time
clear all
close all
clc
set(0,'defaulttextinterpreter', 'tex');

%% Import data
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\filmstrip_selection';
load_fname_str = 'travelling_wave_osc_background_network_16_single-double_band_TW';
%folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\Fig5-self-organisation';
%load_fname_str = 'sample_complex_trajectory_network_34';

load(fullfile(folder, load_fname_str), 'cells_hist', 'positions', 'distances', 'save_consts_struct');
rcell = save_consts_struct.rcell;
a0 = save_consts_struct.a0;
t_out = numel(cells_hist)-1;

% calculate periodicity
[period_ub, t_onset] = periodicity_test_short(cells_hist);
[period, t_onset] = periodicity_test_detailed(cells_hist, t_onset,...
    period_ub);

% SPECIFY NETWORK
network = 34;
appendix = ''; % Note special rule for 33a, 33b
networks_all = [15 19 33 33 34 36];
network_idx = find(network==networks_all, 1);
if strcmp(appendix, 'b')
    network_idx = 4; % special case: 33b
end

% specify other parameters
save_folder = 'H:\My Documents\Multicellular automaton\figures\TW_predictors';

%% Plot wave score (power spectra for 2 genes)
t_start = 0;
t_stop = numel(cells_hist)-1;
[power_spec_dif] = calc_power_spec_dif(t_start, t_stop, cells_hist);

h=figure;

plot(t_start:t_out, power_spec_dif);
%xlim([0 4000]);
ylim([0 1]);
xlabel('t');
ylabel('Wave score (2 genes)');
set(gca, 'FontSize', 20);

fname_str = sprintf('wave_score_2_genes_vs_t_%s', load_fname_str);
fname = fullfile(save_folder, fname_str);
qsave = 1;
save_figure(h, 10, 8, fname, '.pdf', qsave);
%% Plot wave score (power spectra for 4 states)
[wave_dist_vec, mean_dist_vec, background] = determine_wave_dist(t_start,t_stop,cells_hist);

h=figure;
plot(t_start:t_out, wave_dist_vec);
%xlim([0 4000]);
ylim([0 1]);
xlabel('t');
ylabel('Wave score (4 states)');
set(gca, 'FontSize', 20);

fname_str = sprintf('wave_score_4_states_vs_t_%s', load_fname_str);
fname = fullfile(save_folder, fname_str);
qsave = 1;
save_figure(h, 10, 8, fname, '.pdf', qsave);