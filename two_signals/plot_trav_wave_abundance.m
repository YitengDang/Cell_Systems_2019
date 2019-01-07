%% Plots the travelling wave abundances calculated from analytical estimate (Mathematica data)
clear all
%% Load data
folder = 'H:\My Documents\Multicellular automaton\data\two_signals';
fname_str = 'NwaveDensity'; %'TravWaveDensity';
fname = fullfile(folder, fname_str);
xls_data = xlsread(fname);
%% Plot theoretical estimates

idx_sel = 1:size(xls_data, 2);

h = figure;
plot(xls_data(1,idx_sel), 1./xls_data(2,idx_sel), 'bo-');
set(gca, 'YScale', 'log');
xlabel('Grid size');
%ylabel('Abundance');
ylabel('Time');
set(gca, 'FontSize', 24);

%% Plot simulation results
%gz_all = [8:2:20 25]; % from simulations

% New version that loads simulation data can be found in:
% -> analyze_trajectories_trav_wave_vs_variable 