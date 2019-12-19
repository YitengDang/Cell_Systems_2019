% Plots moving averages (also variances and cv)
clear all
close all
clc

%% Load trajectory data
load_folder = 'H:\My Documents\Multicellular automaton\paper_2\figures\originals\Fig5-self-organisation\5A_sample_trajectories_all';
%fname_str = '5B_trajectory_ini_p1_0p4_p2_0p4_v1_no_TW_formed';
%load(fullfile(load_folder, fname_str));

[fname_str, path] = uigetfile(fullfile(load_folder, '\*.mat'), 'Load saved simulation');
load(fullfile(path, fname_str));

% set save folder
%save_folder = 'H:\My Documents\Multicellular automaton\paper_2\figures\data\Fig_5_trajectories';
save_folder = 'H:\My Documents\Multicellular automaton\paper_2\figures\originals\Fig5-self-organisation\5A_sample_trajectories_all';

s = size(cells_hist{1}, 2);
N = size(cells_hist{1}, 1);
tmax = numel(cells_hist)-1;
a0 = save_consts_struct.a0;

%% p(t)
% calculate p(t)
Non = zeros(tmax+1, s);

for i=1:tmax+1
    cells = cells_hist{i};
    for j=1:s
        Non(i,j) = sum(cells(:, j));
    end
end

% Plot p(t)
pt = Non/N;
fig_pos = [1 1 10 8];
t0 = 0;
h = plot_p_vs_t(cells_hist(1:21), t0, fig_pos);
box on

% Save p(t) plot
[~, fname_str_2, ~] = fileparts(fname_str);
fname_str_save = sprintf('%s_p_vs_t_0to20', fname_str_2);
fname_out = fullfile(save_folder, fname_str_save);
qsave = 0;
save_figure(h, 9, 4, fname_out, '.pdf', qsave)
% insets: 9 x 4

%% calculate moving qts p(t)
k = 100; %window length
p_mov_mean = movmean(pt,k,1);
p_mov_var = movvar(pt,k,0,1);
p_mov_cv = p_mov_var./p_mov_mean;
%%
h=figure;
box on
hold on
plot_clrs = [1 0 0;
            0 0 1];
%plot(0:tmax, p_mov_mean);
%plot(0:tmax, p_mov_var);
for i=1:s
    clr = plot_clrs(i, :);
    %p1=plot(0:tmax, p_mov_cv(:, i), 'Color', clr, 'LineWidth', 2);
    p1=plot(0:200, p_mov_cv(1:201, i), 'Color', clr, 'LineWidth', 2);
end
p1.Color(4) = 0.75; %transparency
xlabel('t', 'Interpreter', 'tex');
ylabel('CV(p) = \sigma_p/\mu_p', 'Interpreter', 'tex');
%legend(num2cell(string(1:s)));
set(gca, 'FontSize', 24);
%xlim([0 tmax])
xlim([0 150])
ylim([0 0.02]);

%{
% plot together with pt
yyaxis right
for i=1:s
    clr = plot_clrs(i, :);
    plot(0:tmax, pt(:, i), '--', 'Color', clr);
end
%}
        
[~, fname_str_2, ~] = fileparts(fname_str);
fname_str_save = sprintf('%s_moving_avg_pt_window_%d_t0to150', fname_str_2, k);
fname_out = fullfile(save_folder, fname_str_save);
qsave = 1;
save_figure(h, 9, 4, fname_out, '.pdf', qsave)
%save_figure(h_fig, width, height, path_out, ext, qsave, colored_background)

%% I(t)
% Get data
It = zeros(tmax+1, s);
Theta = zeros(tmax+1, s);
% fixed positions
for i=1:tmax+1
    for j=1:s
        cells = cells_hist{i};
        %[I(i,j), Theta(i,j)] = moranI(cells, a0*dist);
        [It(i,j), Theta(i,j)] = moranI(cells(:, j), a0*distances);
    end
end
%{
% moving cells
for i=1:tmax+1
    for j=1:s
        % calculate new distance each time step
        gz = sqrt(N);
        Lx = 1;
        delx = Lx/gz;
        dely = sqrt(3)/2*delx;
        Ly = dely*gz;
        pos = pos_hist{i};
        dist = calc_dist_periodic(pos(:,1), pos(:,2), Lx, Ly);

        % calculate I, Theta
        %cells = cells_hist{i}{j};
        cells = cells_hist{i};
        %[I(i,j), Theta(i,j)] = moranI(cells, a0*dist);
        [I(i,j), Theta(i,j)] = moranI(cells(:, j), a0*dist);
    end
end
%}

%% I(t) averaging
k = 100; %window length
I_mov_mean = movmean(It,k,1);
I_mov_var = movvar(It,k,0,1);
I_mov_cv = I_mov_var./I_mov_mean;

%% plot I(t)
t0 = 0;
pos_hist = [];
option = 1; % 1: I, 2: Theta
fig_pos = [1 1 9 4];
h = plot_I_vs_t(cells_hist(1:21), t0, a0, distances, option, fig_pos);
%h = plot_I_vs_t(cells_hist(1:51), t0, a0, distances, option, fig_pos);
ylim([0 0.6]);
box on
% plot_I_vs_t(cells_hist, t0, a0, dist, option, fig_pos)

% Save I(t) plot
[~, fname_str_2, ~] = fileparts(fname_str);
fname_str_save = sprintf('%s_I_vs_t_initial_t_0_20', fname_str_2);
fname_out = fullfile(save_folder, fname_str_save);
qsave = 0;
save_figure(h, fig_pos(3), fig_pos(4), fname_out, '.pdf', qsave)
%%
% Plot CV vs t
h=figure;
box on
hold on
plot_clrs = [1 0 0;
            0 0 1];
%plot(0:tmax, p_mov_mean);
%plot(0:tmax, p_mov_var);
for i=1:s
    clr = plot_clrs(i, :);
    %p1 = plot(0:tmax, I_mov_cv(:, i), 'Color', clr, 'LineWidth', 3);
    p1 = plot(0:200, I_mov_cv(1:201, i), 'Color', clr, 'LineWidth', 1.2);
end
p1.Color(4) = 0.75; %transparency
xlabel('t', 'Interpreter', 'tex');
ylabel('CV(I) = \sigma_I/\mu_I', 'Interpreter', 'tex');
%legend(num2cell(string(1:s)));
set(gca, 'FontSize', 24);
%xlim([0 tmax])
xlim([0 150])
ylim([0 0.03]);

[~, fname_str_2, ~] = fileparts(fname_str);
fname_str_save = sprintf('%s_moving_avg_It_window_%d_t0to150', fname_str_2, k);
fname_out = fullfile(save_folder, fname_str_save);
qsave = 1;
save_figure(h, 9, 4, fname_out, '.pdf', qsave)
%save_figure(h_fig, width, height, path_out, ext, qsave, colored_background)