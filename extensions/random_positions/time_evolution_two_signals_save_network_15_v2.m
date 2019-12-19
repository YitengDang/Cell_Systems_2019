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
% * use the randomization algorithm to place the cells on a
% different lattice
% * inner loop over K12 to keep the number of simulations with given
% parameters more or less constant
% * decrease the number of periodicity checks to one every t_check time
% steps
% * shortened loop over simulations by defining a function for each run
% input_ini_state: manually input an initial state (from excel file)
% * run simulations with exact same initial state, only vary extension
% parameters
close all
clear all
rng('shuffle');
%% Simulation parameters
remote = 0;

% variable to loop over
% sigma_D_all = 10.^[-3 -2 -1];
mcsteps_all = [1 10 100 1000];
%noise_all = [0.001 0.01 0.1 1];
%noise_all = [0.01 0.03 0.05 0.1 0.3 0.5 1];
%noise_all = [0.002 0.02 0.2]; %[0.002 0.02 0.2];

% number of simulations to do 
sim_count = 500;

% other settings
network = 15;
networks_all = [15 19 33 34 36];
tmax = 10^4; % max. number of time steps 
% InitiateI = 0; % generate lattice with input I?

% folder to save simulations in
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\randomized lattice';
if remote
    parent_folder = strrep(parent_folder, 'N:\', 'W:\staff-bulk\');
end
subfolder = sprintf('TW_formation_network_%d', network);
save_folder = fullfile(parent_folder, subfolder);

% default file name
sim_ID = 'two_signal_mult';

% number of parameters
num_params = 10; % min( size(Con_wave_sim, 1), 10 ); % number of parameter sets to simulate

%% (2) Load parameters that spontaneously generate TWs from batch simulations 
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2';
fname_str = 'batch_sim_all_topologies_run2_count_TW_analyzed';
if remote
    folder = strrep(folder, 'N:\', 'W:\staff-bulk\');
end

% load K, Con parameters
load(fullfile(folder, fname_str), 'Con_all_by_network', 'K_all_by_network');
network_idx = find(networks_all==network, 1);

Con_wave_sim = Con_all_by_network{network_idx};
K_wave_sim = K_all_by_network{network_idx};

Con_wave_sim = unique(Con_wave_sim, 'rows');
num_psets = size(Con_wave_sim, 1);
K_wave_sim = unique(K_wave_sim(:,:), 'rows');
K_wave_sim = reshape(K_wave_sim, num_psets, 2, 2);

fprintf('Number of parameter sets: %d \n',...
    num_psets);

% load other parameters
fname_str = 'batch_sim_analyzed_data_batch2';
load(fullfile(folder, fname_str), 'N', 'a0', 'hill', 'noise', 'lambda',...
    'InitiateI', 'rcell', 'M_int_all_reduced');

% choose random Con, K values
% Con = Con_wave_sim(1,:);
% K = squeeze(K_wave_sim(1,:,:));

lambda12 = lambda(2);
Coff = [1 1];
Rcell = rcell*a0;
gz = sqrt(N);
M_int = M_int_all_reduced{network};

%InitiateI = 0;
cell_type = zeros(N,1);

mcsteps = 0;
growth_rate = 0;
R_division = 0;
sigma_D = 0;

p0 = Inf;
I0 = Inf;

% Generate lattice
nodisplay = 1;
[positions, distances] = initial_cells_random_markov_periodic(gz, mcsteps,...
    rcell, nodisplay);
%% Calculate # required simulations
% Loop over Con, K values
sim_to_do = zeros(num_params, numel(mcsteps_all));

% folder
%{
folder = save_folder;
if exist(folder, 'dir') ~= 7
    warning('Folder does not exist! ');
    mkdir(folder);
    fprintf('Made new folder %s \n', folder);
end
%}

for idx_param_loop=1:num_params
    Con = Con_wave_sim(idx_param_loop,:);
    K = squeeze(K_wave_sim(idx_param_loop,:,:));

    % folder
    folder = fullfile(save_folder, sprintf('Parameter_set_%d', idx_param_loop));
    if exist(folder, 'dir') ~= 7
        warning('Folder does not exist! ');
        mkdir(folder);
        fprintf('Made new folder %s \n', folder);
    end
        
    for idx_loop=1:numel(mcsteps_all)
        %sigma_D = sigma_D_all(idx_loop);
        mcsteps = mcsteps_all(idx_loop);
        % noise = noise_all(idx_loop);       

%(!!!)  % Filename pattern (!!!)
        pattern = strrep(sprintf('%s_N%d_params_%d_sim_%s_mcsteps_%d_t_out_%s_period_%s',...
            sim_ID, N, idx_param_loop, '(\d+)', mcsteps, '(\d+)', '(\d+|Inf)' ),...
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
        fprintf('N=%d, parameter set %d, mcsteps = %d, sim to do: %d \n',...
            N, idx_param_loop, mcsteps, sim_count-filecount);
        
        sim_to_do(idx_param_loop, idx_loop) = sim_count-filecount;
    end
end
fprintf('Total number of simulations to do: %d \n', sum(sim_to_do(:)) );
%% Then, do the simulations
for idx_param_loop=1:num_params
    %% save folder
    this_save_folder = fullfile(save_folder, sprintf('Parameter_set_%d', idx_param_loop));
    
    % ===================== get original simulations ======================
    % folder for original simulations (without extension)
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_reliability\TW_formation_network_15';
    subfolder = sprintf('Parameter_set_%d', idx_param_loop);
    folder_orig_sim = fullfile(load_folder, subfolder);
    
    % get all files in folder
    listing = dir(folder_orig_sim);
    num_files = numel(listing)-2;
    names_orig_sim = {};
    filecount = 0;
    pattern = strrep(sprintf('%s_N%d_ini_state_rand_params_%d_sim_%s_t_out_%s_period_%s',...
        sim_ID, N, idx_param_loop, '(\d+)', '(\d+)', '(\d+|Inf)' ),...
        '.', 'p');
    for i = 1:num_files
        filename = listing(i+2).name;
        [~, name, ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            [tokens, match] = regexp(filename, pattern, 'tokens', 'match');
            %disp(match);
            if ~isempty(match)
                filecount = filecount + 1;
                file_idx = str2double(tokens{1}{1});
                names_orig_sim{file_idx} = name;
                %names_orig_sim{end+1} = name;
            end
        end
    end
    
    if ~all(~cellfun(@isempty, names_orig_sim)) || numel(names_orig_sim)<500
        error('File names error');
    end
    %%
    % =====================================================================
    for trial=1:sim_count
        %% ============ load initial state =================================
        this_name = names_orig_sim{trial};
        this_fname = fullfile(folder_orig_sim, this_name);
        temp = load(this_fname, 'cells_hist', 'save_consts_struct');
        cells_ini = temp.cells_hist{1};
        K = temp.save_consts_struct.K;
        Con = temp.save_consts_struct.Con;
        %%
        % copy file to new folder
        %{
        new_fname = fullfile(this_save_folder, this_name);
        [status, msg] = copyfile(this_fname, new_fname);
        if ~status
            error(msg);
        end
        %}
        % =================================================================        
        for idx_loop=1:numel(mcsteps_all)
            %sigma_D = sigma_D_all(idx_loop);
            mcsteps = mcsteps_all(idx_loop);
            
            % skip simulation if enough simulations have been done
            filecount = sim_count-sim_to_do(idx_param_loop, idx_loop);
            if trial <= filecount
                continue;
            end
            
            fprintf('Param. set %d, mcsteps %d, trial %d \n', idx_param_loop, mcsteps, trial);
            
            % ----------- simulation ------------------------------------
            display_fig = 0;
            positions = {};
            distances = {};
            fname_str_template = strrep(sprintf('%s_N%d_params_%d_sim_%d_mcsteps_%d',...
                sim_ID, N, idx_param_loop, trial, mcsteps), '.', 'p');
            
            [cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks(...
                N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
                distances, positions, sim_ID, mcsteps, InitiateI, p0, I0, cells_ini,...
                tmax, this_save_folder, fname_str_template, display_fig);
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
    
    if exist(folder, 'dir') ~= 7
        warning('Folder does not exist! ');
        break
    end
    
    % Filename pattern
    % !!!
    %pattern = strrep(sprintf('%s_sigma_D_%.3f_t_out_%s_period_%s-v%s',...
    %	sim_ID, sigma_D, '(\d+)', '(\d+|Inf)', '(\d+)'), '.', 'p');
    pattern = strrep(sprintf('%s_N%d_ini_state_rand_mcsteps_%d_t_out_%s_period_%s',...
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