%% Batch simulations of large number of randomly generated parameter sets
% For 1 fixed topology (specified by M_int)
% Load parameters from earlier set
clear variables
close all
clc
% maxNumCompThreads(6); % limits the number of cores used by MATLAB

%% Settings
%n_pset = 10000; % number of parameter sets to do
nsim = 100; % number of simulations per parameter set
remote = 0; %remote server?

% wave characteristics
num_waves = 1;
wave_type = 1;

% Network
network = 15; % manually set network number (for saving files)
appendix = '';

% get index of network
network_all = [15 19 33 33 34 36];
network_idx = find(network == network_all, 1);
if strcmp(appendix, 'b')
    network_idx = network_idx+1; % exception for network 33 (two wave forms)
end

% Fixed parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];
lambda12  = lambda(2)/lambda(1);
hill = Inf;
noise = 0;
Coff = [1 1];

% Initial conditions
p0 = Inf; % set to Inf = random initial conditions
I0 = Inf; 
tmax = 10000;
InitiateI = 0;

% get pos, dist
mcsteps = 0;
[positions, distances] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Folder for storing simulations
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\TW_formation_with_prop_parameters';

% other
sim_ID = 'two_signal_mult';
cells_ini = []; % []=randomly generate

wave_forms = [3     4     2     1;
     4     3     1     2;
     4     2     1     3;
     3     4     2     1;
     4     2     1     3;
     2     4     3     1];
states_perm = wave_forms(network_idx, :);

%% Get M_int
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2';
fname_str = 'batch_sim_all_topologies_M_int_all_reduced';
load(fullfile(folder, fname_str));
%M_int_all = M_int_all_reduced([15 19 33 33 34 36]);
M_int = M_int_all_reduced{network};

%% Process raw data -> extract parameters for wave propagation (simulations) -> save in separate files
%{

% Get wave forms & permutations (for file naming)
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
if remote
    folder = strrep(folder, 'N:\', 'W:\staff-bulk\');
end
subfolder = 'run2';
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_%s',...
    num_waves, wave_type, subfolder);
load( fullfile(folder, fname_str) );

% get waveforms and networks of possible waves
[x_found, y_found] = find(wave_possible);
num_psets = numel(x_found);
t=table(P(x_found, :), y_found, 'VariableNames', {'Wave_type', 'Network'});
t2=table(x_found, y_found, 'VariableNames', {'wave_idx', 'Network'});

disp('Cell states F, M, B, E');
disp(t);
disp(t2);
for idx_loop=1:6
%idx_loop = 1;

wave_idx = x_found(idx_loop);
network = y_found(idx_loop);
fprintf('wave_idx %d, network %d \n', wave_idx, network);
states_perm = P(wave_idx, :);

load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2';
if remote
    load_folder = strrep(load_folder, 'N:\', 'W:\staff-bulk\');
end
fname_str = sprintf('Trav_wave_predictor_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d_processed',...
       num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4)); 
fname = fullfile(load_folder, fname_str);

load(fname, 'Con_all', 'K_all', 'trav_wave_conds_met');  % load parameter sets and data

%%
% load simulation data
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2_stability_sim';
if remote
    load_folder = strrep(load_folder, 'N:\', 'W:\staff-bulk\');
end
fname_str = sprintf('stability_sim_from_pred_trav_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d',...
       num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));       
fname = fullfile(load_folder, fname_str);
load(fname, 'trav_wave_all', 'trav_wave_all_strict', 'unmodified_all');

% get list of simulations that generated trav. waves
trav_wave_conds_met = (trav_wave_all & unmodified_all);
fprintf('Network %d, #trav. wave conditions = %d \n', network, sum(trav_wave_conds_met) )
    
% filter data on TW simulations
Con_wave_sim = Con_all(trav_wave_conds_met, :);
K_wave_sim = K_all(trav_wave_conds_met, :, :);

%% Save K, Con data in separate file
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
fname_str = sprintf('TW_prop_parameters_simulation_K_Con_network_%d_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
    network, num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
save(fullfile(save_folder, fname_str), 'K_wave_sim', 'Con_wave_sim', 'N', 'a0', 'Coff', 'M_int', 'hill',...
    'noise', 'p0', 'I0', 'rcell', 'Rcell',  'lambda12', 'lambda', ...
    'sim_ID', 'InitiateI', 'tmax', 'mcsteps')

end
%}
%% Load parameter ssets
% ---- Load selected random parameter set ----
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
fname_str = sprintf('TW_prop_parameters_simulation_K_Con_network_%d_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
    network, num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
load(fullfile(load_folder, fname_str), 'K_wave_sim', 'Con_wave_sim');

% ---- Generate random parameter set ----
%{
% bounds
K_b = [1 10^3]; 
Con_b = [1 10^3];

% Latin hypercube
idxK = find(M_int~=0);
idxCon = find(sum(abs(M_int), 1)~=0);

nK = numel(idxK); %sum(sum(abs(M_int))); 
nCon = numel(idxCon); %sum(sum(abs(M_int), 1)>0);

x = lhsdesign(n_pset, nK+nCon);
%}
% ----------------------------------------

% Visualize parameters
%{
figure;
hold on
xlim(K_b);
ylim(Con_b);
%}

%% Loop over parameter sets
for idx_pset=1:size(Con_wave_sim, 1)
    %{
    % from generated data
    thisK = zeros(2);
    thisCon = zeros(1,2);
    thisK(idxK) = (K_b(2) - K_b(1))*x(idx_pset, 1:nK) + K_b(1);
    thisCon(idxCon) = (Con_b(2) - Con_b(1))*x(idx_pset, nK+1:end) + Con_b(1); 
    %}
    
    % from loaded data
    thisK = squeeze(K_wave_sim(idx_pset, :, :));
    thisCon = Con_wave_sim(idx_pset, :);
    
    % Visualize parameters
    % plot(thisK(2,2), thisCon(2), 'bo');
    % plot(thisK(1,2), thisK(2,1), 'bo');
    % plot(thisCon(1), thisCon(2), 'ro');

    % get save folder
    subfolder1 = sprintf('Network_%d_N%d', network, N);
    subfolder2 = sprintf('Param_%d', idx_pset);
    save_folder = fullfile(parent_folder, subfolder1, subfolder2);

    if exist(save_folder, 'dir')~=7
        mkdir(save_folder);
    end

    % save / load parameter set
    fname = fullfile(save_folder, 'parameters.mat');
    if exist(fname, 'file')==2
        try 
            load(fname, 'thisK', 'thisCon');
        catch ME
            if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
                warning('Could not find parameters.mat file for topology %d, pset %d',...
                    network, idx_pset);
                break
            end
        end
    else
        save(fname, 'thisK', 'thisCon');
    end
    
    % count # simulations to do
    %pattern = 'all_topologies_simulate-v(\d+)';
    pattern = 'Simulate_network_(\d+)_params_(\d+)_t_out_(\d+)_period_(\d+|Inf)-v(\d+)';
    [sim_todo, ~] = batch_sim_all_topologies_count_todo(...
        nsim, save_folder, pattern);
    fprintf('Parameter set %d, sims to do: %d \n', idx_pset, sim_todo);
    
    % additional input
    fname_str_template = sprintf('Simulate_network_%d_params_%d', network, idx_pset);
    display_fig = 0;
    
    % simulate trajectories
    %
    for count=1:sim_todo
        %disp(count);
        [cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks(...
            N, a0, Rcell, lambda, hill, noise, M_int, thisK, thisCon, Coff,...
            distances, positions, sim_ID, mcsteps, InitiateI, p0, I0, cells_ini, tmax,...
            save_folder, fname_str_template, display_fig);
    end
    %}
end

%% Test parameter sets
%{
% Load initial conditions
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_snapshots';
%fname_str = 'trav_wave_single_horizontal_inward_bend';
%fname_str = 'trav_wave_single_horizontal_outward_bend'; 
fname_str = 'trav_wave_single_vertical'; 

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

%%
idx = 1352;
Con = Con_wave_sim(idx,:);
K = squeeze(K_wave_sim(idx,:,:));
tmax = 100;
disp_mol = 12;
showI = 0;

cells = cells_in;

cells_hist = {};
cells_hist{end+1} = cells;
    
hin=figure;
plot_handle = reset_cell_figure(hin, positions, rcell);
t = 0;
period = Inf; %default values
t_onset = Inf; 
[cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
        Rcell, Con, Coff, K, lambda, hill, noise);
update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, distances);

% always check within first t_ac time steps
t_ac = 10^2; 
while changed && period==Inf && t<t_ac
    %disp(t);
    pause(0.1);
    t = t+1;
    cells = cellsOut;
    cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
    [period, t_onset] = periodicity_test_short(cells_hist); 
    update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, distances);
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, distances, M_int, a0,...
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