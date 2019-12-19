% Runs multiple simulations of the system with two signalling molecules and saves
% them.
% Parameters can be input by hand, or loaded from a saved simulation file
% with the correct syntax.
% First counts how many simulations to run, based on the number of
% simulation files already present in the specified folder.
% Then, runs the simulations by looping over some given parameter loop, in
% such a way that each parameter set should have the same number of final
% simulations (+-1).

% ----------- Updates ----------------
% v3: use the randomization algorithm to place the cells on a
% different lattice
% v4: inner loop over K12 to keep the number of simulations with given
% parameters more or less constant
% v5: decrease the number of periodicity checks to one every t_check time
% steps
% v6: shortened loop over simulations by defining a function for each run
% input_ini_state: manually input an initial state (from excel file)
close all
clear all
%% Parameters and settings
% TU desktop or remote server?
remote = 0;

% variable to loop over
% sigma_D_all = 10.^[-3 -2 -1];
mcsteps_all = [10 100 1000];

% number of simulations to do 
sim_count = 3;

% other settings
% InitiateI = 0; % generate lattice with input I?
networks_all = [15 19 33 33 34 36]; 
network = 33; 
appendix = 'b'; % Note special rule for 33a, 33b
tmax = 10^4; % max. number of time steps 

network_idx = find(network==networks_all, 1);
if strcmp(appendix, 'b')
    network_idx = 4; % special case: 33b
end

% folder to save simulations in
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\randomized lattice';
if remote
    parent_folder = strrep(parent_folder, 'N:\', 'W:\staff-bulk\');
end

%subfolder = sprintf('K12_%d', K12_all);
subfolder = sprintf('TW_propagation_network_%d%s', network, appendix);
save_folder = fullfile(parent_folder, subfolder);
            
% default file name
sim_ID = 'two_signal_mult';

%% (2) Load parameters from TW stability analysis
% Load simulation parameters
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2_net_parameters_TW_sim';
if remote
    folder = strrep(folder, 'N:\', 'W:\staff-bulk\');
end
fname_str_all = {...
    'Wave_type_1_network_15_states_F3_M4_B2_E1_Con_K_values_waves_sim';
    'Wave_type_1_network_19_states_F4_M3_B1_E2_Con_K_values_waves_sim';
    'Wave_type_1_network_33_states_F3_M4_B2_E1_Con_K_values_waves_sim';
    'Wave_type_1_network_33_states_F4_M2_B1_E3_Con_K_values_waves_sim';
    'Wave_type_1_network_34_states_F4_M2_B1_E3_Con_K_values_waves_sim';
    'Wave_type_1_network_36_states_F2_M4_B3_E1_Con_K_values_waves_sim';
    };
fname_str = fname_str_all{network_idx};
load(fullfile(folder, fname_str), 'N', 'a0', 'hill', 'lambda', 'noise',...
    'rcell', 'Con_wave_sim', 'K_wave_sim');
disp(fname_str);

% Load M_int from another file 
folder2 = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2';
fname_str2 = 'batch_sim_analyzed_data_batch2';
if remote
    folder2 = strrep(folder2, 'N:\', 'W:\staff-bulk\');
end
load(fullfile(folder2, fname_str2), 'M_int_all_reduced');

% choose random Con, K values
%Con = Con_wave_sim(1,:);
%K = squeeze(K_wave_sim(1,:,:));

lambda12 = lambda(2);
Coff = [1 1];
Rcell = rcell*a0;
gz = sqrt(N);
M_int = M_int_all_reduced{network};

InitiateI = 0;
cell_type = zeros(N,1);

growth_rate = 0;
R_division = 0;
sigma_D = 0;
%% Load initial state
signal_count = 2;
%folder = 'D:\Multicellularity\app\data\system_states';
folder = 'H:\My Documents\Multicellular automaton\app\data\system_states';
if remote
    folder = strrep(folder, 'H:\', 'W:\staff-homes\d\yitengdang\');
end

fname = fullfile(folder, 'trav_wave_single_vertical_central_position');
[status, cells_ini, ini_state_fname] = manual_input_state(signal_count, folder, N, fname);

cell_states_all = cell(6, 1);
cell_states_all{1} = [1 0; 1 1; 0 1; 0 0]; % F, M, B, E
cell_states_all{2} = [1 1; 1 0; 0 0; 0 1]; % F, M, B, E
cell_states_all{3} = cell_states_all{1}; % F, M, B, E
cell_states_all{4} = [1 1; 0 1; 0 0; 1 0]; % F, M, B, E
cell_states_all{5} = cell_states_all{4}; % F, M, B, E
cell_states_all{6} = [0 1; 1 1; 1 0; 0 0]; % F, M, B, E
cell_states = cell_states_all{network_idx};

cells_idx00_E = find( cells_ini*[1; 2]==0 );
cells_idx01_F = find( cells_ini*[1; 2]==1 );
cells_idx10_B = find( cells_ini*[1; 2]==2 );
cells_idx11_M = find( cells_ini*[1; 2]==3 );

cells_ini(cells_idx01_F, :) = repmat(cell_states(1,:), gz, 1);
cells_ini(cells_idx11_M, :) = repmat(cell_states(2,:), gz, 1);
cells_ini(cells_idx10_B, :) = repmat(cell_states(3,:), gz, 1);
cells_ini(cells_idx00_E, :) = repmat(cell_states(4,:), N-3*gz, 1);

% generate lattice
nodisplay = 1; 
mcsteps = 0;
[pos_ini, dist_ini] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);
p0 = mean(cells_ini, 1);
I0 = zeros(2,1);
I0(1) = moranI(cells_ini(:,1), a0*dist_ini);
I0(2) = moranI(cells_ini(:,2), a0*dist_ini);

% check initial state
h = figure;
plot_handle = reset_cell_figure(h, pos_ini, rcell);
t=0; disp_mol = 12; showI = 0; 
update_figure_periodic_scatter(plot_handle, cells_ini, t, disp_mol, showI, a0, dist_ini);

%% Calculate # required simulations
% Loop over Con, K values
num_params = min( size(Con_wave_sim, 1), 100 );

sim_to_do = zeros(num_params, numel(mcsteps_all));

% folder
folder = save_folder;
if exist(folder, 'dir') ~= 7
    warning('Folder does not exist! ');
    mkdir(folder);
    fprintf('Made new folder %s \n', folder);
end
        
for idx_param_loop=1:num_params
    Con = Con_wave_sim(idx_param_loop,:);
    K = squeeze(K_wave_sim(idx_param_loop,:,:));
    
    for idx_loop=1:numel(mcsteps_all)
        %sigma_D = sigma_D_all(idx_loop);
        mcsteps = mcsteps_all(idx_loop);
        
%(!!!)  % Filename pattern (!!!)
        pattern = strrep(sprintf('%s_N%d_ini_state_TW_params_%d_mcsteps_%d_t_out_%s_period_%s',...
            sim_ID, N, idx_param_loop, mcsteps, '(\d+)', '(\d+|Inf)' ),...
            '.', 'p');

        listing = dir(folder);
        num_files = numel(listing)-2;
        names = {};
        filecount = 0;
        for i = 1:num_files
            filename = listing(i+2).name;
            % remove extension and do not include txt files
            [~,name,ext] = fileparts(filename);
            if strcmp(ext, '.mat')
                match = regexp(name, pattern, 'match');
                %disp(match);
                if ~isempty(match)
                    filecount = filecount + 1;
                    names{end+1} = name;
                end
            end
        end

        %fprintf('N=%d, sigma_D = %.2f sim to do: %d \n', N, sigma_D, sim_count-filecount);
        fprintf('N=%d, parameter set %d, mcsteps = %d sim to do: %d \n',...
            N, idx_param_loop, mcsteps, sim_count-filecount);
        sim_to_do(idx_param_loop, idx_loop) = sim_count-filecount;
    end

end
fprintf('Total number of simulations to do: %d \n', sum(sim_to_do(:)) );

%% Then, do the simulations
for idx_param_loop=1:num_params
    Con = Con_wave_sim(idx_param_loop,:);
    K = squeeze(K_wave_sim(idx_param_loop,:,:));

for idx_loop=1:numel(mcsteps_all)
    %sigma_D = sigma_D_all(idx_loop);
    mcsteps = mcsteps_all(idx_loop);
    for trial=1:sim_count
        fprintf('Param. set %d, mcsteps %d, trial %d \n', idx_param_loop, mcsteps, trial);

        % skip simulation if enough simulations have been done
        if trial > sim_to_do(idx_param_loop, idx_loop)
            continue;
        end
        % ----------- simulation ------------------------------------
        display_fig = 0;
        positions = {};
        distances = {};
        fname_str_template = strrep(sprintf('%s_N%d_ini_state_TW_params_%d_mcsteps_%d',...
        	sim_ID, N, idx_param_loop, mcsteps), '.', 'p');
        
        [cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks(...
            N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
            distances, positions, sim_ID, mcsteps, InitiateI, p0, I0, cells_ini,...
            tmax, save_folder, fname_str_template, display_fig);
        %}
        %{
        [cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks_moving_cells(...
            N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
            distances, positions, mcsteps, sigma_D, cells_ini, ...
            growth_rate, R_division, sim_ID, tmax, save_folder, display_fig);
        %}
        %--------------------------------------------------------------------------
        %}
    end
end

end
%% (1) Load parameters from saved trajectory
%{
% with parameters saved as structure array 
% load data
data_folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\moving_cells';
file = 'subdomain_oscillation_sigmaD_0_neg_control';
%[file, data_folder] = uigetfile(fullfile(data_folder, '\*.mat'), 'Load saved simulation');
load(fullfile(data_folder, file));

s = save_consts_struct;
N = s.N;
a0 = s.a0;
K = s.K;
Con = s.Con;
Coff = s.Coff;
M_int = s.M_int;
hill = s.hill;
noise = s.noise;
rcell = s.rcell;
lambda12 = s.lambda12;
lambda = [1 lambda12];
mcsteps = str2double(s.mcsteps);

p0 = s.p_ini;
%tmax =  s.tmax;
gz = sqrt(N);
Rcell = rcell*a0;

cell_type = zeros(N,1);

% Initial I
InitiateI = 0;
I0 = [0 0];
s_fields = fieldnames(s);
for i=1:numel(s_fields)
    if strcmp(s_fields{i},'I_ini_str')
        if ~isempty(s.I_ini_str)
            I0 = s.I_ini;
            InitiateI = 1;
        end
    end
end
%

I_ini_str = '';
if InitiateI
    I_ini_str = sprintf('_I_ini_%.2f_%.2f', I0(1), I0(2));
end
%}

%% First, calculate how many simulations are needed 
% without inner loop

%{
sim_to_do = zeros(numel(mcsteps_all));

for idx_loop=1:numel(mcsteps_all)
    %sigma_D = sigma_D_all(idx_loop);
    mcsteps = mcsteps_all(idx_loop);
    
    % subfolder
    folder = save_folder;
    if exist(folder, 'dir') ~= 7
        mkdir(folder);
    end
    %
    if exist(folder, 'dir') ~= 7
        warning('Folder does not exist! ');
        break
    end
    
    % Filename pattern
    % !!!
    %pattern = strrep(sprintf('%s_sigma_D_%.3f_t_out_%s_period_%s-v%s',...
    %	sim_ID, sigma_D, '(\d+)', '(\d+|Inf)', '(\d+)'), '.', 'p');
    pattern = strrep(sprintf('%s_N%d_ini_state_TW_mcsteps_%d_t_out_%s_period_%s',...
    	sim_ID, N, mcsteps, '(\d+)', '(\d+|Inf)' ), '.', 'p');
    
    listing = dir(folder);
    num_files = numel(listing)-2;
    names = {};
    filecount = 0;
    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            match = regexp(name, pattern, 'match');
            %disp(match);
            if ~isempty(match)
                filecount = filecount + 1;
                names{end+1} = name;
            end
        end
    end

    %fprintf('N=%d, sigma_D = %.2f sim to do: %d \n', N, sigma_D, sim_count-filecount);
    fprintf('N=%d, mcsteps = %d sim to do: %d \n', N, mcsteps, sim_count-filecount);
    sim_to_do(idx_loop) = sim_count-filecount;
end
%}