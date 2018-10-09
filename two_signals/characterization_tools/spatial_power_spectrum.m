% Computes the spatial power spectrum of a macroscopic trajectory
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');

%% Load trajectory
folder = 'K:\bn\hy\Shared\Yiteng\Multicellularity\videos\selected\data';
%folder = 'K:\bn\hy\Shared\Yiteng\Multicellularity\videos\synchronization';
load_fname = 'plane_wave_formation_period_15';
load(fullfile(folder, load_fname), 'cells_hist', 'positions', 'distances', 'save_consts_struct');

% Process loaded data
t_out = numel(cells_hist)-1;
N = size(cells_hist{1}, 1);
gz = sqrt(N);
a0 = save_consts_struct.a0;
rcell = save_consts_struct.rcell;

% Default save folder 
%save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\fourier_analysis';
save_folder = 'K:\bn\hy\Shared\Yiteng\Multicellularity\videos\selected\spatial power spectra';
%% Save a lattice snapshot 
% For comparison with spectral density, save a lattice snapshot
time = 400; % time of the simulation to look at
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

%% 2D Fourier transform
cells_norm = (cells - mean(cells, 1)); % normalise cells

cells_fft1 = fft2(reshape(cells_norm(:,1), gz, gz));
%cells_fft1(1,1) = 0;
cells_fft2 = fft2(reshape(cells_norm(:,2), gz, gz));
%cells_fft2(1,1) = 0;

power1 = abs(fftshift(cells_fft1)).^2/N;
power2 = abs(fftshift(cells_fft2)).^2/N;

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
qsave = 1;
fname_str = sprintf('%s_spatial_spectrum_t%d', load_fname, time);
fname = fullfile(save_folder, fname_str);
colored_background = 0;
save_figure(h, 0, 0, fname, '.pdf', qsave, colored_background);
%% Take long time averages
t1 = 0;
t2 = t_out;

cells_fft1 = zeros(gz, gz);
for t=t1:t2
    cells = cells_hist{t+1};
    cells_norm = (cells - mean(cells, 1)); % normalise cells
    
    cells_fft1 = cells_fft1 + fft2(reshape(cells_norm(:,1), gz, gz));
    cells_fft2 = cells_fft1 + fft2(reshape(cells_norm(:,2), gz, gz));
end
cells_fft1 = fftshift(cells_fft1);
cells_fft2 = fftshift(cells_fft2);
power_lt1 = abs(cells_fft1/(t2-t1+1)).^2/N; % check normalization
power_lt2 = abs(cells_fft2/(t2-t1+1)).^2/N;

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
qsave = 1;
fname_str = sprintf('%s_spatial_spectrum_avg_t%d_to_%d', load_fname, t1, t2);
fname = fullfile(save_folder, fname_str);
colored_background = 0;
save_figure(h, 0, 0, fname, '.pdf', qsave, colored_background);