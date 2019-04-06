clear all
close all
set(0, 'defaulttextinterpreter', 'tex');

%% Parameters
% remote destination (Webdrive)?
remote = 0;

% parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];

%Con = [18 16];
%Coff = [1 1];
%M_int = [0 1; -1 1];
%K = [0 9; 11 6];

hill = Inf;
noise = 0;

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Specify wave type
num_waves = 1;
wave_type = 1;
fname_str_all = {'trav_wave_single_vertical',...
    'trav_wave_single_horizontal_inward_bend',...
    'trav_wave_single_horizontal_outward_bend'};
fname_str = fname_str_all{wave_type};

%% Load overview data
% Load list of found networks and wave forms
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
%%
for idx_loop=1
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    states_perm = P(wave_idx, :);
    
    % load predictor data to get parameter sets
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2';
    if remote
        load_folder = strrep(load_folder, 'N:\', 'W:\staff-bulk\');
    end
    fname_str = sprintf('Trav_wave_predictor_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d_processed',...
           num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4)); 
    fname = fullfile(load_folder, fname_str);
    
    %load(fname);
    load(fname, 'Con_all', 'K_all');
end


%% Analyze data

Con_wave_sim = Con_all(trav_wave_conds_met, :);
K_wave_sim = K_all(trav_wave_conds_met, :, :);
    
P_data = [Con_all]
