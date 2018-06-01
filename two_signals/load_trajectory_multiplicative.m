%% Loads and analyzes a trajectory for the 2 signals, multiplicative system
clear all 
close all

%data_folder = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution';
%fig_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\time_evolution';
%data_folder = 'D:\Multicellularity\data\two_signals\time_evolution';
data_folder = 'D:\Multicellularity\app\Multicellularity-2.1\data\time_evolution';
fig_folder = 'D:\Multicellularity\figures\two_signals\time_evolution';

[file, path] = uigetfile(fullfile(data_folder, '\*.mat'), 'Load saved simulation');
load(fullfile(path, file));

N = save_consts_struct.N;
a0 = save_consts_struct.a0;
gz = sqrt(N);
[dist, pos] = init_dist_hex(gz,gz);
cell_type = zeros(N, 1);
%% replay trajectory
hin = figure();
disp_mol = 1;
for i=1:numel(cells_hist)
    cells = cells_hist{i+1};
    update_cell_figure_continuum(hin, pos, cells, cell_type, i, disp_mol);
end
%%
max_period = periodicity_test_short(cells_hist); % ignore displayed output 
if max_period~=Inf
    [period, t_onset] = periodicity_test_v2(cells_hist);
else 
    period = Inf;
end

if period==Inf
    fprintf('No periodic dynamics found \n');
%else
%    fprintf('Found period %d, starting at t=%d \n',...
%        period, t_onset);
end

%% Plot p vs t
t1 = 1;
t2 = numel(cells_hist)-1;
fig_pos = [1 1 8 6];
msg = plot_p_vs_t(cells_hist(t1:t2), t1, fig_pos);
disp(msg);

% save
qsave = 0;
if qsave
    h = figure(1);
    [~, fname_str] = fileparts(file);
    label = sprintf('_p_vs_t_t%d_%d', t1, t2);
    path_out = fullfile(fig_folder, ...
        strcat(fname_str, label));
    ext = '.pdf';
    save_figure(h, 10, 8, path_out, ext)
end
%% Plot I vs t
t1 = 5000;
t2 = 5130;
fig_pos = [1 1 8 6];
msg = plot_I_vs_t(cells_hist(t1:t2), t1, a0, dist, 1, fig_pos);
disp(msg);

% save
qsave = 1;
if qsave
    h = figure(2);
    [~, fname_str] = fileparts(file);
    label = sprintf('_I_vs_t_t%d_%d', t1, t2);
    path_out = fullfile(fig_folder, strcat(fname_str, label));
    ext = '.pdf';
    save_figure(h, 10, 8, path_out, ext)
end