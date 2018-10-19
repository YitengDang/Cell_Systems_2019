% Tests whether a travelling wave can propagate by performing explicit
% simulations
clear all
close all
set(0,'defaulttextinterpreter', 'latex')
%% Parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];

Con = [18 16];
Coff = [1 1];
M_int = [0 1; -1 1];
K = [0 9; 11 6];

hill = Inf;
noise = 0;

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% calculate fN
idx = gz + round(gz/2); % pick cell not at corner of grid
dist_vec = a0*dist(idx,:);
r = dist_vec(dist_vec>0);
fN = zeros(2,1);
fN(1) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
fN(2) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));

% save folder
%save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_stability';
%fname_str_default = strrep(sprintf('N%d_a0_%.1f_rcell_%.1f_lambda12_%.1f_M_int_%d_%d_%d_%d_Con_%d_%d',...
%    N, a0, rcell, lambda(2), M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2),...
%    Con(1), Con(2)), '.', 'p');

% Load initial conditions
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_snapshots';
%fname_str = 'trav_wave_single_horizontal_outward_bend'; 
fname_str = 'trav_wave_single_horizontal_outward_bend'; 
%fname_str = 'trav_wave_single_vertical'; 

fname = fullfile(load_folder, fname_str);
cells_load = cell(2,1);
cells_load{1} = xlsread(fname, 'Sheet1');
cells_load{2} = xlsread(fname, 'Sheet2');
if all(size(cells_load{1})==[N 2])
    cells_in = cells_load{1};
elseif all(size(cells_load{1})==[gz gz]) && all(size(cells_load{2})==[gz gz])
    cells_in(:, 1) = reshape(cells_load{1}, N, 1);
    cells_in(:, 2) = reshape(cells_load{2}, N, 1);
else
    disp('Wrong input format');
end

% get cell state population of initial state
cells_idx = cells_in*[1; 2];
n_in = histcounts(cells_idx, -0.5:3.5);

%% Run single simulation (test)
%{
tmax = 100;
disp_mol = 12;
showI = 0;

cells = cells_in;

[trav_wave, cellsOut] = trav_wave_sim_tester(cells_in, pos, dist, rcell, a0, M_int, Con,...
    Coff, K, lambda, hill, noise, tmax);
%}
%{
cells_hist = {};
cells_hist{end+1} = cells;
    
hin=figure;
plot_handle = reset_cell_figure(hin, pos, rcell);
t = 0;
period = Inf; %default values
t_onset = Inf; 
[cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
        Rcell, Con, Coff, K, lambda, hill, noise);
update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);

% always check within first t_ac time steps
t_ac = 10^2; 
while changed && period==Inf && t<t_ac
    %disp(t);
    pause(0.1);
    t = t+1;
    cells = cellsOut;
    cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
    [period, t_onset] = periodicity_test_short(cells_hist); 
    update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
        Rcell, Con, Coff, K, lambda, hill, noise);
end

t_out = t; %default t_out
% if periodicity found, refine check to find period
if period<Inf && t>t_ac
    [period, t_onset] = periodicity_test_detailed(cells_hist, t_check, period);
    t_out = t_onset + period; 
end

if changed && t==tmax
    tmax_string = '_tmax_reached';
else
    tmax_string = '';
end
fprintf('Final: t_out = %d, period %d \n', t_out, period);

[trav_wave] = travelling_wave_test(cells_hist, a0, period, t_out);
fprintf('Travelling wave? %d (1=Yes, 0=No) \n', trav_wave);
%}
%% Check for a set of parameters
%
tmax = 100; % evolve system for tmax steps
K_12_all = 1:30;
K_21_all = 1:30;
K_22_all = 1:10;

trav_wave_all = zeros(numel(K_12_all), numel(K_21_all), numel(K_22_all)); % is it still a trav wave?
unmodified_all = zeros(numel(K_12_all), numel(K_21_all), numel(K_22_all)); % is the travelling wave unmodified (on population level)
for k=1:numel(K_22_all)
    for i=1:numel(K_12_all)
        for j=1:numel(K_21_all)
            thisK = K;
            thisK(1,2) = K_12_all(i);
            thisK(2,1) = K_21_all(j);
            thisK(2,2) = K_22_all(k);
            fprintf('K = [0 %d; %d %d] \n', thisK(1,2), thisK(2,1), thisK(2,2));
            
            [trav_wave, cells_out] = trav_wave_sim_tester(cells_in, pos, dist, rcell, a0, M_int, Con,...
                Coff, thisK, lambda, hill, noise, tmax);
            
            % check whether population has changed
            cells_out_idx = cells_out*[1; 2];
            n_out = histcounts(cells_out_idx, -0.5:3.5);
            fprintf('Unmodified? %d \n', all(n_out==n_in));
            
            % save results
            trav_wave_all(i,j,k) = trav_wave; % whether it is still a travelling wave
            unmodified_all(i,j,k) = all(n_out==n_in); % whether the population of the wave is unmodified
        end        
    end
end
close all
%}

%% Save data
%{
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\data';
fname_str_data = sprintf('%s_stability_sim_K12_K21_K22_range', fname_str);  
fname = fullfile(save_folder, fname_str_data);
save(fname, 'K_12_all', 'K_21_all', 'K_22_all', 'trav_wave_all', 'unmodified_all');
%}

%% Load data
%
load_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\data';
%fname_str_data = sprintf('%s_stability_sim_%s_K12_K21_K22_range', fname_str, fname_str_default);  
fname_str_data = sprintf('%s_stability_sim_K12_K21_K22_range', fname_str);  
fname = fullfile(load_folder, fname_str_data);
load(fname, 'K_12_all', 'K_21_all', 'K_22_all', 'trav_wave_all', 'unmodified_all');
%}
%% Plot 2D results as heatmap
for K22_idx=1:10
    h = figure;
    hold on
    [K_21_mesh, K_12_mesh] = meshgrid(K_21_all, K_12_all);
    %temp1 = squeeze(trav_wave_cond_met(1, :,:,K22_idx));
    idx1 = (trav_wave_all(:, :, K22_idx)==1);
    scatter(K_12_mesh(idx1), K_21_mesh(idx1), 100, 'o');
    idx2 = (unmodified_all(:, :, K22_idx)==1);
    scatter(K_12_mesh(idx2), K_21_mesh(idx2), 100, 'x');

    legend({'trav. wave', 'unmodified'}, 'Location', 'ne');
    xlabel('$K^{(12)}$');
    ylabel('$K^{(21)}$');
    title(sprintf('$K^{(22)} = %.1f $', K_22_all(K22_idx) ));
    set(gca, 'FontSize', 20);
    set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);

    xlim([K_12_all(1)-1 K_12_all(end)+1]);
    ylim([K_21_all(1)-1 K_21_all(end)+1]);

    qsave = 0;    
    fname_str_default = strrep(sprintf('N%d_a0_%.1f_rcell_%.1f_lambda12_%.1f_Con_%d_%d_K22_%d',...
        N, a0, rcell, lambda(2), Con(1), Con(2), K22_idx), '.', 'p');
    fname_str_save = sprintf('%s_stability_sim_%s_vs_K12_K21', fname_str, fname_str_default);    

    fname = fullfile(save_folder, fname_str_save);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
    %close all
end
%% 3D scatter plot
%
[K_21_mesh, K_12_mesh, K_22_mesh] = meshgrid(K_21_all, K_12_all, K_22_all);

h = figure;
hold on
grid on
%
idx1 = squeeze(trav_wave_all(:,:,:)==1);
%s1 = scatter3(K_12_mesh(idx1), K_21_mesh(idx1), K_22_mesh(idx1), 100, 'o', 'filled');
%s1.MarkerFaceAlpha = 0.8;
idx2 = (unmodified_all(:,:,:)==1);
%s2=scatter3(K_12_mesh(idx2), K_21_mesh(idx2), K_22_mesh(idx2), 100, 'd', 'filled');
%s2.MarkerFaceAlpha = 0.8;
idx12 = idx1 & idx2;
s12 = scatter3(K_12_mesh(idx12), K_21_mesh(idx12), K_22_mesh(idx12), 100, 'o', 'filled');
s12.MarkerFaceAlpha = 0.8;

xlim([K_12_all(1)-1 K_12_all(end)+1]);
ylim([K_21_all(1)-1 K_21_all(end)+1]);
zlim([K_22_all(1)-1 K_22_all(end)+1]);
az = -45; el = 30;
view(az,el);
xlabel('$$K^{(12)}$$');
ylabel('$$K^{(21)}$$');
zlabel('$$K^{(22)}$$');
set(gca, 'FontSize', 20);
%}
qsave = 1;    
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability';
fname_str_save = sprintf('%s_stability_sim_3D_scatter', fname_str);
fname = fullfile(save_folder, fname_str_save);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%% functions
%[trav_wave] = trav_wave_sim_tester(cells_in, pos, dist, rcell, a0, M_int, Con,...
%    Coff, K, lambda, hill, noise, tmax);

function [trav_wave, cellsOut] = trav_wave_sim_tester(cells_in, pos, dist, rcell, a0, M_int, Con,...
    Coff, K, lambda, hill, noise, tmax)
    disp_mol = 12;
    showI = 0;
    Rcell = rcell*a0;
    cells = cells_in;
    cells_hist = {};
    cells_hist{end+1} = cells;

    %hin=figure;
    %plot_handle = reset_cell_figure(hin, pos, rcell);
    t = 0;
    period = Inf; %default values
    %t_onset = Inf; 
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);
    %update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);

    % always check within first t_ac time steps
    while changed && period==Inf && t<tmax
        %disp(t);
        pause(0.1);
        t = t+1;
        cells = cellsOut;
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        [period, ~] = periodicity_test_short(cells_hist); 
        %update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);
        [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise);
    end

    t_out = t;
    fprintf('Final: t_out = %d, period %d \n', t_out, period);
    
    % check whether produced pattern is still a travelling wave
    [trav_wave] = travelling_wave_test(cells_hist, a0, period, t_out);
    fprintf('Travelling wave? %d (1=Yes, 0=No) \n', trav_wave);
end