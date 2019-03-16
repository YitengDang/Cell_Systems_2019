%% Loads and analyzes a trajectory for the 2 signals, multiplicative system
clear all
close all
%%
%data_folder = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution';
%fig_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\time_evolution';
%data_folder = 'D:\Multicellularity\data\two_signals\time_evolution';
%data_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_p0_K12_9';
data_folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity paper 2\movies';
save_folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity paper 2\movies';

[file, path] = uigetfile(fullfile(data_folder, '\*.mat'), 'Load saved simulation');
load(fullfile(path, file));

N = save_consts_struct.N;
a0 = save_consts_struct.a0;
gz = sqrt(N);
[dist, pos] = init_dist_hex(gz,gz);
cell_type = zeros(N, 1);
rcell = save_consts_struct.rcell;

%% replay trajectory, save "normal" movie
%
savemovie = 1;
if savemovie
    [~, fname_str] = fileparts(file);
    fname_out = fullfile(save_folder, strcat(fname_str, '_simulation')); 
    format = 'Motion JPEG AVI'; %'MPEG-4'; % <- .mp4
    myVideo = VideoWriter(fname_out, format); %, 'Uncompressed AVI');
    myVideo.FrameRate = 10;  % Default 30
    myVideo.Quality = 50; % Default 75
    open(myVideo);
end

h = figure;
disp_mol = 12;
t_out = numel(cells_hist)-1;
%[h_cells, h_borders]  = reset_cell_figure_minimal(h, positions, rcell);
[h_cells, h_borders] = reset_cell_figure(h, positions, rcell);
for i=0:t_out
    pause(0.01);
    cells = cells_hist{i+1};
    %update_cell_figure_continuum(hin, pos, cells, cell_type, i, disp_mol);
    %update_figure_periodic_scatter(plot_handle, cells, time, disp_mol, showI, a0, distances)
    update_cell_figure_external(h_cells, h_borders, cells, i, disp_mol, positions);    
    
    if savemovie
        frame = getframe(h);
        writeVideo(myVideo, frame);
    end
end
if savemovie
    close(myVideo);
end
%}
%% Check periodicity (does not always work well)
[period, t_onset] = periodicity_test_short(cells_hist); 
t_check_init = 1;
if period<Inf
    [period, t_onset] = periodicity_test_detailed(cells_hist, t_check_init, period);
    t_out = t_onset + period;
    fprintf('Final: t_out = %d, period %d \n', t_out, period);
end

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
folder = 'H:\My Documents\Multicellular automaton\paper_2\figures\data\Fig_5_trajectories\all_5B';
[~, fname_str_pre, ~] = fileparts(file);
fname_str = sprintf('%s_p1_p2_dynamics', fname_str_pre); 
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
p01 = scatter( Non(end,1)/N, Non(end,2)/N, 100, 'rs', 'filled');
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);
xlim([0 1])
ylim([0 1])
xlabel('p^{(1)}', 'Interpreter', 'tex');
ylabel('p^{(2)}', 'Interpreter', 'tex');
legend([p00 p01], {'Initial state', 'Final state'}, 'AutoUpdate','off');
set(gca, 'FontSize', 24);
for t_final=0:t2
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
pause(1);
% close video
close(myVideo);