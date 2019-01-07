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
%% Simulation parameters
remote = 0;

% variable to loop over
% sigma_D_all = 10.^[-3 -2 -1];
% mcsteps_all = [0 10 100 1000];
noise_all = [0 0.01 0.05 0.1 0.5 1];
%noise_all = [0.01 0.03 0.05 0.1 0.3 0.5 1];
%noise_all = [0 0.1 1];

% number of simulations to do 
sim_count = 100;

% other settings
network = 15;
networks_all = [15 19 33 34 36];
tmax = 10^4; % max. number of time steps 
% InitiateI = 0; % generate lattice with input I?

% folder to save simulations in
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\\trav_wave_with_noise';
if remote
    parent_folder = strrep(parent_folder, 'N:\', 'W:\staff-bulk\');
end
subfolder = sprintf('TW_formation_network_%d_fixed_parameter_set', network);
save_folder = fullfile(parent_folder, subfolder);
            
% default file name
sim_ID = 'two_signal_mult';

%% Load fixed parameters
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\sweep K12 new lattice\N225';
fname_str = 'two_signal_mult_N225_initiateI0_randpos_mcsteps0_K12_9_t_out_78_period_15-v1';
load(fullfile(load_folder, fname_str), 'save_consts_struct');

N = save_consts_struct.N;
a0 = save_consts_struct.a0;
hill = save_consts_struct.hill;
noise = save_consts_struct.noise;
lambda12 = save_consts_struct.lambda12;
rcell = save_consts_struct.rcell;
M_int = save_consts_struct.M_int; 
Con = save_consts_struct.Con;
K = save_consts_struct.K;

%lambda12 = lambda(2);
lambda = [1 lambda12];
Coff = [1 1];
Rcell = rcell*a0;
gz = sqrt(N);

InitiateI = 0;
cell_type = zeros(N,1);

mcsteps = 0;
growth_rate = 0;
R_division = 0;
sigma_D = 0;

p0 = Inf;
I0 = Inf;
%% Load initial state
%{
signal_count = 2;
%folder = 'D:\Multicellularity\app\data\system_states';
folder = 'H:\My Documents\Multicellular automaton\app\data\system_states';
fname = fullfile(folder, 'trav_wave_single_vertical_central_position');
[status, cells_ini, ini_state_fname] = manual_input_state(signal_count, folder, N, fname);

nodisplay = 1; 
mcsteps = 0;
[~, dist_ini] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);
p0 = mean(cells_ini, 1);
I0 = zeros(2,1);
I0(1) = moranI(cells_ini(:,1), a0*dist_ini);
I0(2) = moranI(cells_ini(:,2), a0*dist_ini);
%}
%% Calculate # required simulations
sim_to_do = zeros(numel(noise_all), 1);

% folder
folder = save_folder;
if exist(folder, 'dir') ~= 7
    warning('Folder does not exist! ');
    mkdir(folder);
    fprintf('Made new folder %s \n', folder);
end

for idx_loop=1:numel(noise_all)
    %sigma_D = sigma_D_all(idx_loop);
    % mcsteps = mcsteps_all(idx_loop);
    noise = noise_all(idx_loop);

    %(!!!)  % Filename pattern (!!!)
    pattern = strrep(sprintf('%s_N%d_ini_state_rand_fixed_params_noise_%.3f_t_out_%s_period_%s',...
        sim_ID, N, noise, '(\d+)', '(\d+|Inf)' ),...
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
    %fprintf('N=%d, parameter set %d, mcsteps = %d sim to do: %d \n',...
    %    N, idx_param_loop, mcsteps, sim_count-filecount);
    fprintf('N=%d, noise = %.3f, sim to do: %d \n',...
        N, noise, sim_count-filecount);

    sim_to_do(idx_loop) = sim_count-filecount;
end

fprintf('Total number of simulations to do: %d \n', sum(sim_to_do(:)) );

%% Then, do the simulations
for idx_loop=1:numel(noise_all)
    %sigma_D = sigma_D_all(idx_loop);
    %mcsteps = mcsteps_all(idx_loop);
    noise = noise_all(idx_loop);
    for trial=1:sim_count
        %fprintf('Param. set %d, mcsteps %d, trial %d \n', idx_param_loop, mcsteps, trial);
        fprintf('Noise %.3f, trial %d \n', noise, trial);
        
        % skip simulation if enough simulations have been done
        if trial > sim_to_do(idx_loop)
            continue;
        end
        % ----------- simulation ------------------------------------
        display_fig = 0;
        positions = {};
        distances = {};
        fname_str_template = strrep(sprintf('%s_N%d_ini_state_rand_fixed_params_noise_%.3f',...
        	sim_ID, N, noise), '.', 'p');
        cells_ini = [];
        
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