% Computes the spatial power spectrum of the snapshot of the system
% Get a snapshot of the system by loading a full trajectory and looking at
% the system at a particular time. Decompose the lattice of the snapshot
% into Fourier modes and calculate the spatial power spectrum. Take a
% long-time average of the power spectra over the entire simulation.
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');

%% Load trajectory
folder = 'K:\bn\hy\Shared\Yiteng\Multicellularity\videos\selected\data';
%folder = 'K:\bn\hy\Shared\Yiteng\Multicellularity\videos\synchronization';
load_fname = 'plane_wave_formation_period_15';
load_fname = 'Chaos_in_growing_domain';
load(fullfile(folder, load_fname), 'cells_hist', 'positions', 'distances', 'save_consts_struct');

% Process loaded data
t_out = numel(cells_hist)-1;
N = size(cells_hist{1}, 1);
gz = sqrt(N);
a0 = save_consts_struct.a0;
rcell = save_consts_struct.rcell;

% Default save folder for figures
%save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\fourier_analysis';
save_folder = 'K:\bn\hy\Shared\Yiteng\Multicellularity\videos\selected\spatial power spectra';
%% Load and snapshot of the lattice
% For comparison with spectral density, save a lattice snapshot
% First, specify the time of the simulation to look at
time = 400; 
cells = cells_hist{time+1};

% Plot snapshot
cell_type = zeros(N,1);
disp_mol = 12;
showI = 0;
h = figure;
plot_handle = reset_cell_figure(h, positions, rcell);
update_figure_periodic_scatter(plot_handle, cells, time, disp_mol, showI, a0, distances)
%update_cell_figure_continuum(h, positions, cells, cell_type, time, disp_mol)
set(gcf, 'color', 'w');

% Save snapshot
qsave = 0;
fname_str = sprintf('%s_snapshot_t%d', load_fname, time);
fname = fullfile(save_folder, fname_str);
colored_background = 1;
save_figure(h, 0, 0, fname, '.pdf', qsave, colored_background);

%% Plot the spatial power spectrum (2D discrete Fourier transform)
% 2D Fourier transform, calculate power spectrum

% normalise cells
cells_norm = zeros(size(cells));
cells_norm(:, 1) = (cells(:, 1) - mean(cells(:, 1)))./std(cells(:, 1))/gz; 
cells_norm(:, 2) = (cells(:, 2) - mean(cells(:, 2)))./std(cells(:, 2))/gz; 

cells_fft1 = fft2(reshape(cells_norm(:,1), gz, gz));
cells_fft2 = fft2(reshape(cells_norm(:,2), gz, gz));

power1 = abs(fftshift(cells_fft1)).^2/(N-1);
power2 = abs(fftshift(cells_fft2)).^2/(N-1);

% Plot power spectrum
h = figure;
x = 0:gz-1;
x = fftshift(0:gz-1);
idx = find(x==0);
x(1:idx-1) = x(1:idx-1) - gz;
imagesc(x/gz, x/gz, power2)
xlabel('$q_x$');
ylabel('$q_y$');
title(sprintf('t=%d', time));
c=colorbar;
ylabel(c, '$$S(q_x, q_y)$$', 'Interpreter', 'latex');
set(gca, 'YDir', 'normal', 'FontSize', 20)
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);

% Save figure
qsave = 0;
fname_str = sprintf('%s_spatial_spectrum_t%d', load_fname, time);
fname = fullfile(save_folder, fname_str);
colored_background = 0;
save_figure(h, 0, 0, fname, '.pdf', qsave, colored_background);
%% Check normalisation snapshot
%{
cells_norm = zeros(size(cells));
cells_norm(:, 1) = (cells(:, 1) - mean(cells(:, 1)))./std(cells(:, 1))/gz; % normalise cells
cells_norm(:, 2) = (cells(:, 2) - mean(cells(:, 2)))./std(cells(:, 2))/gz; % normalise cells
mean(cells_norm)
var(cells_norm)

cells_fft1 = fft2(reshape(cells_norm(:,1), gz, gz));
%cells_fft2 = fft2(reshape(cells_norm(:,2), gz, gz));
mean(cells_fft1(:))
sum(sum(abs(cells_fft1).^2))/N
%}
%% Take long time averages
t1 = 0;
t2 = t_out;

cells_fft1 = zeros(gz, gz);
cells_fft2 = zeros(gz, gz);
cells_norm_all = {};
cells_fft1_all = {};
cells_fft2_all = {};
for t=t1:t2
    cells = cells_hist{t+1};
    
    % normalise cells
    cells_norm = zeros(size(cells));
    cells_norm(:, 1) = (cells(:, 1) - mean(cells(:, 1)))./std(cells(:, 1))/gz; 
    cells_norm(:, 2) = (cells(:, 2) - mean(cells(:, 2)))./std(cells(:, 2))/gz; 

    cells_fft1 = cells_fft1 + fft2(reshape(cells_norm(:,1), gz, gz));
    cells_fft2 = cells_fft1 + fft2(reshape(cells_norm(:,2), gz, gz));
    
    cells_norm_all{end+1} = cells_norm;
    cells_fft1_all{end+1} = fft2(reshape(cells_norm(:,1), gz, gz));
    cells_fft2_all{end+1} = fft2(reshape(cells_norm(:,2), gz, gz));
end
cells_fft1 = fftshift(cells_fft1);
cells_fft2 = fftshift(cells_fft2);
power_lt1 = abs(cells_fft1/(t2-t1+1)).^2/(N-1); % check normalization
power_lt2 = abs(cells_fft2/(t2-t1+1)).^2/(N-1);

h=figure;
x = fftshift(0:gz-1);
idx = find(x==0);
x(1:idx-1) = x(1:idx-1) - gz;
imagesc(x, x, power_lt1)
xlabel('$q_x$');
ylabel('$q_y$');
c=colorbar;
ylabel(c, '$$S(q_x, q_y)$$', 'Interpreter', 'latex');
set(gca, 'YDir', 'normal', 'FontSize', 20)
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);
title('Long-time average');

% Save figure
qsave = 0;
fname_str = sprintf('%s_spatial_spectrum_avg_t%d_to_%d', load_fname, t1, t2);
fname = fullfile(save_folder, fname_str);
colored_background = 0;
save_figure(h, 0, 0, fname, '.pdf', qsave, colored_background);

%% Check normalisation over time
nsteps = numel(cells_fft1_all);
fft_sum2 = zeros(nsteps, 2); 
cells_sum2 = zeros(nsteps, 2);
for i=1:nsteps
    fft_sum2(i, 1) = sum(abs(cells_fft1_all{i}(:)).^2);
    fft_sum2(i, 2) = sum(abs(cells_fft2_all{i}(:)).^2);
    cells_sum2(i,1) = sum(cells_norm_all{i}(:,1).^2);
    cells_sum2(i,2) = sum(cells_norm_all{i}(:,2).^2);  
end

figure;
hold on
plot(1:nsteps, fft_sum2/(N-1));  % spectrum normalisation
plot(1:nsteps, cells_sum2*N/(N-1)); % cell state normalisation
ylim([0.8 1.2]);
