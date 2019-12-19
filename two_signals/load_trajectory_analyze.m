%% Loads and analyzes a trajectory for the 2 signals, multiplicative system
clear all 
close all
%%
%data_folder = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution';
%fig_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\time_evolution';
%data_folder = 'D:\Multicellularity\data\two_signals\time_evolution';
%data_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_p0_K12_9';
%data_folder = 'H:\My Documents\Multicellular automaton\paper_2\figures\originals\Fig5-self-organisation\5A_sample_trajectories_all';
%data_folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution';
data_folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity paper 2\movies\other_for_presentations';
save_folder = data_folder;

[file, path] = uigetfile(fullfile(data_folder, '\*.mat'), 'Load saved simulation');
load(fullfile(path, file));

N = save_consts_struct.N;
a0 = save_consts_struct.a0;
gz = sqrt(N);
[dist, pos] = init_dist_hex(gz,gz);
cell_type = zeros(N, 1);
rcell = save_consts_struct.rcell;
if ~exist('pos_hist', 'var')
    pos_hist = {};
end
%% Save trajectory
[~, fname_raw, ~] = fileparts(file);%'trajectory_ini_p1_0p4_p2_0p4_v1_TW_formed';
fname = fullfile(save_folder, fname_raw);
%save(fname, 'cells_hist', 'positions', 'distances', 'save_consts_struct')

%% replay trajectory
h = figure;
disp_mol = 12;

[h_cells, h_borders]  = reset_cell_figure_minimal(h, positions, rcell);  
for i=0:numel(cells_hist)-1
    pause(0.1);
    cells = cells_hist{i+1};
    %update_cell_figure_continuum(hin, pos, cells, cell_type, i, disp_mol);
    %update_figure_periodic_scatter(plot_handle, cells, time, disp_mol, showI, a0, distances)
    update_cell_figure_external(h_cells, h_borders, cells, i, disp_mol, positions);    
    
    %{
    frame = getframe(gcf);
    img =  frame2im(frame);
    [img,cmap] = rgb2ind(img,256);
    if i == 0
        imwrite(img,cmap,'animation.gif','gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(img,cmap,'animation.gif','gif','WriteMode','append','DelayTime',1);
    end
    %}
end

%% Save movie (.gif)
qsave = 1;
disp_mol = 12;
fname_out = fullfile(save_folder, file);
%fname_out = 'H:\My Documents\Multicellular automaton\temp\temp_jpeg';
frame_rate = 5;
format = 'Motion JPEG AVI';
t_ini = 20;
t_out = 30;
if qsave
    save_movie(cells_hist, rcell, positions, pos_hist, disp_mol, fname_out,...
        frame_rate, format, t_ini, t_out)
    %save_movie_gif(cells_hist, rcell, positions, pos_hist, disp_mol, fname_out,...
    %    frame_rate);
    %save_movie_gif(cells_hist, rcell, positions, pos_hist, disp_mol, fname_out,...
    %    frame_rate, t_ini, t_out)
end

%% Check periodicity (does not always work well)
[period, t_onset] = periodicity_test_short(cells_hist); 
t_check_init = 1;
if period<Inf
    [period, t_onset] = periodicity_test_detailed(cells_hist, t_check_init, period);
    t_out = t_onset + period;
end
fprintf('Final: t_out = %d, period %d \n', t_out, period);
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
    path_out = fullfile(save_folder, ...
        strcat(fname_str, label));
    ext = '.pdf';
    save_figure(h, 10, 8, path_out, ext)
end
%% Plot I vs t
t1 = 1;
t2 = numel(cells_hist)-1;
fig_pos = [1 1 8 6];
msg = plot_I_vs_t(cells_hist(t1:t2), t1, a0, dist, 1, fig_pos);
box on
disp(msg);

% save
qsave = 0;
if qsave
    h = figure(2);
    [~, fname_str] = fileparts(file);
    label = sprintf('_I_vs_t_t%d_%d', t1, t2);
    path_out = fullfile(save_folder, strcat(fname_str, label));
    ext = '.pdf';
    save_figure(h, 10, 8, path_out, ext, qsave)
end

%% Plot p1, p2 vs t
% calculate Non(t)
t1 = 0;
%t2 = 300;
t2 = numel(cells_hist)-1;

Non = zeros(t2-t1+1, 2);
for i=t1+1:t2+1
    cells = cells_hist{i};
    for j=1:2
        Non(i,j) = sum(cells(:, j));
    end
end

fig_pos = [1 1 8 6];

h = figure();
cla(h, 'reset');
hold on
% plot selection
t_TW = 1687;
dt = 100;
t_idx = (0:dt) + 1;
%t_idx = (tt-dt):tt +1;
plot( Non(t_idx,1)/N, Non(t_idx,2)/N, 'b-' );
scatter( Non(t_idx(1), 1)/N, Non(t_idx(1), 2)/N, 100, 'gs', 'filled' );
scatter( Non(t_idx(end),1)/N, Non(t_idx(end),2)/N, 100, 'ro', 'filled');
title(sprintf('First %d time steps', dt));
% plot everything
%plot( Non(:,1)/N, Non(:,2)/N, 'b-' );
%scatter( Non(1, 1)/N, Non(1, 2)/N, 100, 'gs', 'filled' );
%scatter( Non(end,1)/N, Non(end,2)/N, 100, 'ro', 'filled');
xlabel('p^{(1)}', 'Interpreter', 'tex');
ylabel('p^{(2)}', 'Interpreter', 'tex');
set(gca, 'FontSize', 24);
xlim([0 1])
ylim([0 1])
box on

qsave = 0;
save_folder = 'H:\My Documents\Multicellular automaton\paper_2\figures\data\Fig_5_trajectories';
%fname_raw = 'trajectory_ini_p1_0p5_p2_0p5_v1_TW_formed';
save_fname_str = sprintf('%s_t_out_%d_plot_p1_p2_vs_t_first_%d_ts', fname_raw, t2, dt);
save_fname = fullfile(save_folder, save_fname_str);
width = 10; 
height = 8;
save_figure(h, width, height, save_fname, '.pdf', qsave)

%% Plot trajectory in p1, p2 space as movie
% Data
t1 = 0;
%t2 = 300;
t2 = numel(cells_hist)-1;
Non = zeros(t2-t1+1, 2);
for i=t1+1:t2+1
    cells = cells_hist{i};
    for j=1:2
        Non(i,j) = sum(cells(:, j));
    end
end

% filename 
folder = 'H:\My Documents\Multicellular automaton\paper_2\figures\data\Fig_5_trajectories';
[~, fname_str_pre, ~] = fileparts(file);
fname_str = sprintf('%s_video', fname_str_pre); 
fname_out = fullfile(folder, fname_str);

% Video settings
frame_rate = 8; % frames/second
format = 'MPEG-4'; %'Motion JPEG AVI'; %movie format 
% 'Motion JPEG AVI' <- default, works best
% 'Uncompressed AVI' <- high quality(?), large file
% 'MPEG-4' <- .mp4
% 'Archival' <- unknown ext
% 'Motion JPEG 2000' <- unknown ext

% open video
myVideo = VideoWriter(fname_out, format); %, 'Uncompressed AVI');
myVideo.FrameRate = frame_rate;  % Default 30
myVideo.Quality = 50; % Default 75
open(myVideo);

% plot frames
t_TW = t2; % formation time TW
h=figure;
box on
%clf(h, 'reset');
p00 = scatter( Non(1, 1)/N, Non(1, 2)/N, 100, 'gs', 'filled');
hold on
p01 = scatter( Non(end,1)/N, Non(end,2)/N, 100, 'ro', 'filled');

set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);
xlim([0 1])
ylim([0 1])
xlabel('p^{(1)}', 'Interpreter', 'tex');
ylabel('p^{(2)}', 'Interpreter', 'tex');
legend([p00 p01], {'Initial state', 'Final state'}, 'AutoUpdate','off');
set(gca, 'FontSize', 24);
for t_final=0:1700
    p1 = plot( Non(1:(t_final+1),1)/N, Non(1:(t_final+1),2)/N, 'b-' );
    p1.Color(4) = 0.6;
    p2 = scatter( Non(t_final+1,1)/N, Non(t_final+1,2)/N, 100, 'ko', 'filled');
    uistack(p00,'top'); % move plots to top
    uistack(p01,'top');
    
    title(sprintf('Time = %d', t_final));
    %if t_final>=t_TW
    %    title(sprintf('Time = %d (TW formed)', t_final));
    %else
    %    title(sprintf('Time = %d', t_final));
    %end
    
    frame = getframe(h);
    writeVideo(myVideo, frame);
    pause(0.001);
    delete(p1);
    delete(p2);
end
hold off 

% close video
close(myVideo);