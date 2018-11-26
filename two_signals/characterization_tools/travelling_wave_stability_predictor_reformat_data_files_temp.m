clear all
close all
set(0,'defaulttextinterpreter', 'latex')

%% Parameters
% Number of parameter sets to do
n_pset = 10^5;

% Manual input
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];
Coff = [1 1];

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% specify wave type and characteristics
wave_types_str = {'straight', 'inward bend', 'outward bend'};
wave_type = 3;
num_waves = 1; % number of waves
bandwidth = 1; % width of band of cells of one type

folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
subfolder = 'run2';
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_%s',...
    num_waves, wave_type, subfolder);
load( fullfile(folder, fname_str) );
%}

%% display found waves
%
[x_found, y_found] = find(wave_possible);
t=table(P(x_found, :), y_found, 'VariableNames', {'Wave_type', 'Network'});
t2=table(x_found, y_found, 'VariableNames', {'wave_idx', 'Network'});

disp('Cell states F, M, B, E');
disp(t);
disp(t2);

%%
for idx_loop=1:numel(x_found)
    %wave_idx = 2;
    %network = 19;
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    
    states_perm = P(wave_idx, :);
    %
    % load data
    load_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
    %subfolder = 'run3_vary_a0_lambda12';
    fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
            num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    load( fullfile(load_folder, subfolder, fname_str) );
    
    trav_wave_conds_met = squeeze(trav_wave_conds_met(network, :));
    K_all = squeeze(K_all(network, :, :, :));
    Con_all = squeeze(Con_all(network, :, :));
    
    % save data
    save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2';
    fname_str = sprintf('Trav_wave_predictor_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d_processed',...
        num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    fname = fullfile(save_folder, fname_str);
    save(fname, 'gz', 'a0', 'rcell', 'lambda',...
        'trav_wave_conds_met', 'Con_all', 'K_all');
end