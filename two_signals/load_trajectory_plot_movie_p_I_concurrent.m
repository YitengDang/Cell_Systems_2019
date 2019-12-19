%% Loads and analyzes a trajectory for the 2 signals, multiplicative system
% Plot the dynamic movie of the macroscopic variables
clear all 
close all
%% Load simulation
data_folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity paper 2\movies\selected for paper\data files';
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

%% Plot trajectory with other graphs together
h = figure;
disp_mol = 12;

[h_cells, h_borders]  = reset_cell_figure_minimal(h, positions, rcell);  
for i=0:numel(cells_hist)-1
    pause(0.1);
    cells = cells_hist{i+1};
    %update_cell_figure_continuum(hin, pos, cells, cell_type, i, disp_mol);
    %update_figure_periodic_scatter(plot_handle, cells, time, disp_mol, showI, a0, distances)
    update_cell_figure_external(h_cells, h_borders, cells, i, disp_mol, positions);    
    
    %
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
qsave = 0;
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

%% Calculate p and I as function of t
tmax = length(cells_hist)-1;
s = 2;
% p(t)
Non = zeros(tmax+1, s);
for i=1:tmax+1
    cells = cells_hist{i};
    for j=1:s
        Non(i,j) = sum(cells(:, j));
    end
end
% I(t)
I_t = zeros(tmax+1, s);
for i=1:tmax+1
    for j=1:s
        cells = cells_hist{i};
        [I_t(i,j), Theta_t(i,j)] = moranI(cells(:, j), a0*dist);
    end
end

%% Generate large plot with subplots

% Video information
folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity paper 2\movies\Fig_5_self-organization';
[~, fname_str_pre, ~] = fileparts(file);
fname_str = sprintf('%s_sim_p1_p2_combined_video_fps5', fname_str_pre); 
fname_out = fullfile(folder, fname_str);

% Video settings
frame_rate = 5; % frames/second
format = 'MPEG-4'; %movie format 
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

% Plot figure, save frames
h = figure;
% loop over times
for t=1:tmax+1
    % ---------- Plot cells ----------
    cells = cells_hist{t};
    h1 = subplot(6, 2, 3:2:9);
    axis off
    set(gcf, 'Color', [0.8 0.8 0.8]);
    %[h_cells, h_borders]  = reset_cell_figure_minimal(h1, positions, rcell);  
    hold on
    % colours
    c_all = ones(N, 3);
    clr_k = zeros(N, 3); % black boundaries
    %markers = {'o', 's'};
    Rcell_px = 5.5;
    h_cells = scatter(positions(:,1), positions(:,2), (2*Rcell_px)^2, c_all, 'filled', 'o');

    % Adjust cell colors
    clrs = 1-cells;
    c_all = zeros(size(cells, 1), 3); 
    c_all(:, 1) = clrs(:, 2);               % signal 2 absent -> Turn on red channel
    c_all(:, 2) = clrs(:, 1).*clrs(:, 2);   % greenness nonlinear function of redness and blueness
    c_all(:, 3) = clrs(:, 1);               % signal 1 absent -> Turn on blue channel
    set(h_cells, 'cdata', c_all);
    title(sprintf('Time = %d', t), 'FontSize', 16);

    % ---------- plot p(t) ---------- 
    h2 = subplot(6, 2, [2 4 6]);
    box on
    hold on
    xlim([0 tmax]);
    ylim([0 1]);
    %xlabel('t', 'Interpreter', 'tex');
    set(gca, 'XTick', []);
    ylabel('p^{(i)}', 'Interpreter', 'tex');
    set(gca, 'FontSize', 16);
    p1 = plot(0:t-1, Non(1:t, 1)/N, 'r', 'LineWidth', 1.5);
    p2 = plot(0:t-1, Non(1:t, 2)/N, 'b', 'LineWidth', 1.5);
    p1.Color(4) = 0.5;
    p2.Color(4) = 0.5;
    legend({'1', '2'});
    % ---------- plot I(t) ---------- 
    h3 = subplot(6, 2, [8 10 12]);
    box on
    hold on
    set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);
    xlim([0 tmax]);
    ylim([-0.2 1]);
    xlabel('time', 'Interpreter', 'tex');
    ylabel('I^{(i)}', 'Interpreter', 'tex');
    set(gca, 'FontSize', 16);
    p1 = plot(0:t-1, I_t(1:t, 1), 'r', 'LineWidth', 1.5);
    p2 = plot(0:t-1, I_t(1:t, 2), 'b', 'LineWidth', 1.5);
    p1.Color(4) = 0.5;
    p2.Color(4) = 0.5;
    legend({'1', '2'});
    % ---------------------------------
    set(h, 'Units', 'inches', 'Position', [1 1 16 9]);
    frame = getframe(h);
    writeVideo(myVideo, frame);
    pause(0.1);
end

% Close video
close(myVideo);

%% Plot p(t) as video
% filename for saving
folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity paper 2\movies\Fig_5_self-organization';
[~, fname_str_pre, ~] = fileparts(file);
fname_str = sprintf('%s_p1_p2_vs_t_video', fname_str_pre); 
fname_out = fullfile(folder, fname_str);

% Video settings
frame_rate = 10; % frames/second
format = 'Motion JPEG AVI'; %movie format 
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

% set up figure
h = figure;
box on
hold on
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);
xlim([0 tmax]);
ylim([0 1]);
xlabel('t', 'Interpreter', 'tex');
ylabel('p^{(i)}', 'Interpreter', 'tex');
set(gca, 'FontSize', 24);
% plot frames
for t=1:tmax+1
    p1 = plot(0:t-1, Non(1:t, 1)/N, 'r', 'LineWidth', 1.5);
    p2 = plot(0:t-1, Non(1:t, 2)/N, 'b', 'LineWidth', 1.5);
    p1.Color(4) = 0.5;
    p2.Color(4) = 0.5;
    
    title(sprintf('Time = %d', t-1));
    legend({'1', '2'});
    frame = getframe(h);
    writeVideo(myVideo, frame);
    pause(0.01);
    if t<tmax+1
        delete(p1);
        delete(p2);
    end
end

% close video
close(myVideo);

%% Plot I(t) as video
% filename for saving
folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity paper 2\movies\Fig_5_self-organization';
[~, fname_str_pre, ~] = fileparts(file);
fname_str = sprintf('%s_I1_I2_vs_t_video', fname_str_pre); 
fname_out = fullfile(folder, fname_str);

% Video settings
frame_rate = 10; % frames/second
format = 'Motion JPEG AVI'; %movie format 
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

% set up figure
h = figure;
box on
hold on
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);
xlim([0 tmax]);
ylim([-0.2 1]);
xlabel('t', 'Interpreter', 'tex');
ylabel('I^{(i)}', 'Interpreter', 'tex');
set(gca, 'FontSize', 24);
% plot frames
for t=1:tmax+1
    p1 = plot(0:t-1, I_t(1:t, 1), 'r', 'LineWidth', 1.5);
    p2 = plot(0:t-1, I_t(1:t, 2), 'b', 'LineWidth', 1.5);
    p1.Color(4) = 0.5;
    p2.Color(4) = 0.5;
    
    title(sprintf('Time = %d', t-1));
    legend({'1', '2'});
    frame = getframe(h);
    writeVideo(myVideo, frame);
    pause(0.01);
    if t<tmax+1
        delete(p1);
        delete(p2);
    end
end

% close video
close(myVideo);