% Save a snapshot of a loaded trajectory
clear all
close all
clc

%% Load trajectory
%folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\Fig_2_Dropbox';
%folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\patterns\Network 33';
%folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\travelling waves';
%folder = 'H:\My Documents\Multicellular automaton\temp';
%folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\filmstrip_selection';
%folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\patterns\Network 12';
%folder = 'H:\My Documents\Multicellular automaton\paper_2\figures\originals\Fig6-extensions\wave_examples';
%folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity paper 2\movies\Fig_6_model_extension_patterns';
folder = 'H:\My Documents\Multicellular automaton\paper_2\figures\originals\Fig5-self-organisation\5A_samples trajectory chosen';
%fname_str = 'travelling_pulse_horizontal_M_int1_1_-1_0_t_out_243_period_30-network_19';
%fname_str = 'spiral_single_period_172_M_int_1_-1_1_0_network_15';
%fname_str = 'sample_complex_trajectory';
%fname_str = 'two_signal_mult_N225_ini_state_TW_params_5_mcsteps_400_t_out_15_period_15-v1';
%fname_str = 'Fig_6_C_TW_formation_hill_10_ini_state_rand';
fname_str = 'sample_complex_trajectory_network_34';

load(fullfile(folder, fname_str), 'cells_hist', 'positions', 'distances', 'save_consts_struct');
%load(fullfile(folder, fname_str), 'cells_hist', 'positions', 'distances', 'save_consts_struct', 'positions_all');
rcell = save_consts_struct.rcell;
a0 = save_consts_struct.a0;

%% Replay trajectory
h = figure;
disp_mol = 12;
[h_cells, h_borders]  = reset_cell_figure_minimal(h, positions, rcell);  
for i=20 %0:numel(cells_hist)-1
    cells = cells_hist{i+1};
    %update_cell_figure_continuum(hin, pos, cells, cell_type, i, disp_mol);
    %update_figure_periodic_scatter(plot_handle, cells, time, disp_mol, showI, a0, distances)
    update_cell_figure_external(h_cells, h_borders, cells, i, disp_mol, positions);    
    
    pause(0.5);
end

%% Plot snapshot(s)
for time=20 %numel(cells_hist)-1
    % parameters
    %time = 27;
    disp_mol = 12;
    showI = 0;

    % plot
    cells = cells_hist{time+1};
    %positions = positions_all{time+1};
    h = figure;
    [h_cells, h_borders]  = reset_cell_figure_minimal(h, positions, rcell);
    %update_figure_periodic_scatter(plot_handle, cells, time, disp_mol, showI, a0, distances)
    update_cell_figure_external(h_cells, h_borders, cells, time, disp_mol, positions);
    %% Save snapshot
    %save_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\Fig_2_Dropbox';
    %save_folder = 'H:\My Documents\Multicellular automaton\latex\13_pattern_analysis\figures';
    %save_folder = 'H:\My Documents\Multicellular automaton\temp';
    %save_folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\FigS3_class_III_oscillating';
    %save_folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\FigS6_types_of_TWs';
    %save_folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\FigS18_static_patterns';
    %save_folder = 'H:\My Documents\Multicellular automaton\paper_2\figures\originals\Fig6-extensions\wave_examples';
    save_folder = fullfile(folder, 'snapshots');
    
    save_fname_str = sprintf('%s_snapshot_t_%d', fname_str, time);
    save_fname = fullfile(save_folder, save_fname_str);

    width = 0; % set to 0 to keep current width/height
    height = 0;
    saveq = 1;
    colored_background = 1;
    save_figure(h, width, height, save_fname, '.pdf', saveq, colored_background)
end

%% Plot p(t)
t0 = 0;
fig_pos = [1 1 6.5 5];
msg = plot_p_vs_t(cells_hist, t0, fig_pos);
box on
h = gcf;
width = 0; % set to 0 to keep current width/height
height = 0;
colored_background = 0;

save_folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\Fig5-sample complex trajectory';
save_fname_str = sprintf('sample_complex_trajectory_p_vs_t_v2');
save_fname = fullfile(save_folder, save_fname_str);

saveq = 0;
save_figure(h, width, height, save_fname, '.pdf', saveq, colored_background)
%% Plot I(t)
t0 = 0;
option = 1;
fig_pos = [1 1 6.5 5];
msg = plot_I_vs_t(cells_hist, t0, a0, distances, option, fig_pos);
ylim([-1 1]);
box on
h = gcf;
width = 0; % set to 0 to keep current width/height
height = 0;
colored_background = 0;

save_folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\Fig5-sample complex trajectory';
save_fname_str = sprintf('sample_complex_trajectory_I_vs_t_v2');
save_fname = fullfile(save_folder, save_fname_str);

saveq = 0;
save_figure(h, width, height, save_fname, '.pdf', saveq, colored_background)