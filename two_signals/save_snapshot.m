% Save a snapshot of a loaded trajectory
clear all
close all
clc

%% Load trajectory
%folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\Fig_2_Dropbox';
%folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\patterns\Network 33';
%folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\travelling waves';
folder = 'H:\My Documents\Multicellular automaton\temp';
%fname_str = 'travelling_pulse_horizontal_M_int1_1_-1_0_t_out_243_period_30-network_19';
%fname_str = 'spiral_single_period_172_M_int_1_-1_1_0_network_15';
fname_str = 'temp_trajectory';
load(fullfile(folder, fname_str), 'cells_hist', 'positions', 'distances', 'save_consts_struct');
rcell = save_consts_struct.rcell;
a0 = save_consts_struct.a0;

%% Plot snapshot
% parameters
time = 0;
disp_mol = 12;
showI = 0;

% plot
cells = cells_hist{time+1};
h = figure;
plot_handle = reset_cell_figure(h, positions, rcell);
update_figure_periodic_scatter(plot_handle, cells, time, disp_mol, showI, a0, distances)

%% Save snapshot
%save_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\Fig_2_Dropbox';
%save_folder = 'H:\My Documents\Multicellular automaton\latex\13_pattern_analysis\figures';
save_folder = 'H:\My Documents\Multicellular automaton\temp';
save_fname_str = sprintf('%s_snapshot_t_%d', fname_str, time);
save_fname = fullfile(save_folder, save_fname_str);

width = 0; % set to 0 to keep current width/height
height = 0;
saveq = 1;
colored_background = 1;
save_figure(h, width, height, save_fname, '.pdf', saveq, colored_background)