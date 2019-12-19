% Save a snapshot of a loaded trajectory
clear all
close all
clc

%% Load trajectory
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\one_signal_finite_Hill\sample_trajectories';
fname_str = 'N121_a0_5p60_K8_Con16_hill2p00_t172_xmeanf_0p69_prec_8_uniform-v1_partial_synchronization';
load(fullfile(load_folder, fname_str), 'cells_hist', 'positions', 'distances', 'save_consts_struct');
%load(fullfile(folder, fname_str), 'cells_hist', 'positions', 'distances', 'save_consts_struct', 'positions_all');
rcell = save_consts_struct.rcell;
a0 = save_consts_struct.a0;

%% Replay trajectory
h = figure;
disp_mol = 1;
[h_cells, h_borders]  = reset_cell_figure_minimal(h, positions, rcell);  
for i=0:numel(cells_hist)-1
    cells = cells_hist{i+1};
    %update_cell_figure_continuum(hin, pos, cells, cell_type, i, disp_mol);
    %update_figure_periodic_scatter(plot_handle, cells, time, disp_mol, showI, a0, distances)
    update_cell_figure_external(h_cells, h_borders, cells, i, disp_mol, positions);    
    
    pause(0.5);
end

%% Plot snapshot(s)

for time=[0 10 20 60 172] %[0:5:20] %numel(cells_hist)-1
    % parameters
    %time = 27;
    disp_mol = 1;
    showI = 0;

    % plot
    cells = cells_hist{time+1};
    %positions = positions_all{time+1};
    h = figure;
    [h_cells, h_borders]  = reset_cell_figure_minimal(h, positions, rcell);
    %update_figure_periodic_scatter(plot_handle, cells, time, disp_mol, showI, a0, distances)
    update_cell_figure_external_temp(h_cells, h_borders, cells, time, disp_mol, positions);
    %% Save snapshot
    save_folder = fullfile(load_folder, 'snapshots');
    
    save_fname_str = sprintf('%s_snapshot_t_%d_GFP', fname_str, time);
    save_fname = fullfile(save_folder, save_fname_str);

    width = 0; % set to 0 to keep current width/height
    height = 0;
    saveq = 1;
    colored_background = 1;
    save_figure(h, width, height, save_fname, '.pdf', saveq, colored_background)
end
%% Color bar (GFP)
rows = 100; cols = 500;
h=figure;
% Plot 1D color map
image=zeros(rows, cols, 3); %initialize
for row=1:rows
    for col=1:cols
        c = (col-1)/(cols-1);
        image(row, col, :) = [0 c 0];
    end
end
imshow(image)
axis on
xlabel('Gene expression');
set(gca, 'XTick', linspace(1, cols, 6), 'XTickLabel', {'0',...
    '0.2', '0.4', '0.6', '0.8', '1'} );
set(gca, 'YTick', [], 'FontSize', 20);    
set(h, 'Position', [100 100 600 200]);

saveq = 0;
colored_background = 0;
save_fname = fullfile(save_folder, 'Colorbar_GFP');
save_figure(h, 8, 3, save_fname, '.pdf', saveq, colored_background)
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