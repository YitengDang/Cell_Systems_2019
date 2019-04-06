% Predicts whether a travelling wave can propagate according to nearest
% neighbour interactions
% Apply computational search over all topologies to try to find parameters
% that allow for propagation of any kind of travelling wave
% rerun: rerun types 2/3 to include the conditions for plane wave
% propagation
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
%M_int = [0 1; -1 1]; % network 15 reversed
%{
M_int = [0 1; -1 1]; % network 15 reversed
M_int = [1 -1; 1 0]; % network 15
M_int = [1 1; -1 0]; % network 19
M_int = [1 -1; 1 1]; % network 33
M_int = [-1 -1; 1 1]; % network 34
M_int = [-1 1; -1 1]; % network 36
%}
%Con = [18 16];
%K = [0 9; 11 4];
%K = [4 11; 9 0];

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Obtain from simulation
%{
folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected';
subfolder = 'patterns\Network 33';
fname_str = 'tripple_wave_diagonally';
fname = fullfile(folder, subfolder, fname_str);
load(fname, 'save_consts_struct');

Con = save_consts_struct.Con;
K = save_consts_struct.K;
N = save_consts_struct.N;
gz = sqrt(N);
a0 = save_consts_struct.a0;
rcell = save_consts_struct.rcell;
Rcell = rcell*a0;
lambda = save_consts_struct.lambda;
M_int = save_consts_struct.M_int;

%}

% specify wave type and characteristics
wave_types_str = {'straight', 'inward bend', 'outward bend'};
wave_type = 3;
num_waves = 1; % number of waves
bandwidth = 1; % width of band of cells of one type

% order of states: F, M, B, E
%states_perm = [3 4 2 1]; % network 15 
%states_perm = [2 4 3 1]; % network 15 reversed
%{
states_perm = [2 4 3 1]; % network 15 reversed
states_perm = [3 4 2 1]; % network 15 
states_perm = [4 3 1 2]; % network 19
states_perm = [3 4 2 1]; % network 33
states_perm = [4 2 1 3]; % network 33/34
states_perm = [2 4 3 1]; % network 36
default_states = [0 0; 0 1; 1 0; 1 1];
%}

% save folder
%subfolder = 'run2'; %run3_vary_a0_lambda12 run2b_n_pset_10e6
subfolder = 'run2_rerun'; 
%save_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2';
save_folder = fullfile('N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general', subfolder);
%% Load analyzed data
%
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_%s',...
    num_waves, wave_type, 'run2');
load( fullfile(load_folder, fname_str) );
%}

% display found waves
[x_found, y_found] = find(wave_possible);
t=table(P(x_found, :), y_found, 'VariableNames', {'Wave_type', 'Network'});
t2=table(x_found, y_found, 'VariableNames', {'wave_idx', 'Network'});

disp('Cell states F, M, B, E');
disp(t);
disp(t2);

% Get a list of all interaction matrices for the found networks
[M_int_found] = get_found_M_int(y_found);

%% Re-test all parameter sets for which we obtained positive results for wave_type=2 or 3
for idx_loop=1:numel(x_found)
    %% Load NNMFA data
    %idx_loop = 1;
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    states_perm = P(wave_idx, :);
    
    % load data
    subfolder = 'run2';
    %fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
    %	num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    fname_str = sprintf('Trav_wave_predictor_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d_processed',...
    	num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4)); 
    fname = fullfile(load_folder, subfolder, fname_str);    
    load( fname );
    
    %% Rerun test for wave_type=1 only
    %{
    % Get relevant parameters
    Con_wave = squeeze(Con_all(network, squeeze(trav_wave_conds_met(network, :)), :));
    K_wave = squeeze(K_all(network, squeeze(trav_wave_conds_met(network, :)), :, :));
    num_wave_sims = size(Con_wave, 1);
    
    trav_wave_conds_met_both_types = zeros(num_wave_sims, 1);
    M_int = M_int_found{idx_loop};
    
    for idx1=1:num_wave_sims
        thisK = squeeze(K_wave(idx1, :, :));
        thisCon = Con_wave(idx1, :);
        thisa0 = a0;
        thislambda = lambda;
        %thisa0 = 10*x(idx1, nK+nCon+1);
        %thislambda = [1 2*x(idx1, nK+nCon+2);];

        this_wave_type = 1;
        trav_wave_conds_met_both_types(idx1) = travelling_wave_stability_predictor_general_func(...
            gz, thisa0, dist, rcell, thislambda, M_int, thisCon, thisK,...
            this_wave_type, states_perm, num_waves, bandwidth);        
    end
    
    idx_both = find(trav_wave_conds_met_both_types);
    Con_wave_new = Con_wave(idx_both, :);
    K_wave_new = K_wave(idx_both, :, :);
    
    % fname
    fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d',...
        num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    fname = fullfile(save_folder, fname_str);
    save(fname, 'gz', 'a0', 'rcell', 'lambda', 'M_int',...
        'Con_wave', 'K_wave', 'trav_wave_conds_met_both_types',...
        'Con_wave_new', 'K_wave_new');
    %}
    %% Rerun test for all wave types
    num_wave_sims = size(Con_all, 1);
    trav_wave_conds_met_both_types = zeros(num_wave_sims, 1);
    M_int = M_int_found{idx_loop};
    
    for idx1=1:num_wave_sims
        thisK = squeeze(K_all(idx1, :, :));
        thisCon = Con_all(idx1, :);
        thisa0 = a0;
        thislambda = lambda;
        %thisa0 = 10*x(idx1, nK+nCon+1);
        %thislambda = [1 2*x(idx1, nK+nCon+2);];

        this_wave_type = 1;
        trav_wave_conds_met_both_types(idx1) = travelling_wave_stability_predictor_general_func(...
            gz, thisa0, dist, rcell, thislambda, M_int, thisCon, thisK,...
            this_wave_type, states_perm, num_waves, bandwidth);        
    end
    
    idx_both = find(trav_wave_conds_met_both_types);
    Con_wave_new = Con_all(idx_both, :);
    K_wave_new = K_all(idx_both, :, :);
    
    % fname
    fname_str = sprintf('trav_wave_conditions_rerun_all_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d',...
        num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    fname = fullfile(save_folder, fname_str);
    save(fname, 'gz', 'a0', 'rcell', 'lambda', 'M_int',...
        'Con_all', 'K_all', 'trav_wave_conds_met_both_types',...
        'Con_wave_new', 'K_wave_new');

end
%% Functions
function [M_int_found] = get_found_M_int(y_found)
    % Get a list of all interaction matrices for the found networks
    M = [0 1 -1]; % index to interaction
    M_int_all_reduced = {};
    done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
    for k=1:3^4
        [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], k);
        gM = [i22 i21; i12 i11];
        M_int = [M(i11) M(i12); M(i21) M(i22)];
        if done(i11,i12,i21,i22)
            continue
        elseif k==1
            continue
        else
            M_int_all_reduced{end+1} = M_int;
            done(i11,i12,i21,i22) = 1;
            done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
        end
    end

    M_int_found = cell(numel(y_found), 1);
    for i=1:numel(y_found)
        M_int_found{i} = M_int_all_reduced{y_found(i)};
        %disp(M_int_found{i})
    end
end