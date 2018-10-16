% v3: use the randomization algorithm to place the cells on a
% different lattice
% v4: inner loop over K12 to keep the number of simulations with given
% parameters more or less constant
% v5: decrease the number of periodicity checks to one every t_check time
% steps
close all
clear all
maxNumCompThreads(3);
%warning off

%% Simulation parameters
max_trials = 100;

% loop 
%p_all = 0:0.1:1; % loop variable (to be specified below)
%noise_all = [0 0.01 0.05 0.1 0.5];
hill_all = [Inf 100 10 5 2];

% other settings
% p0=[0.5 0.5];
InitiateI = 0;
tmax = 10^4;

% folder to save simulations in
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_Hill\K12_9';

%% System parameters
%{
% lattice parameters
%gz = 15;
%N = gz^2;
gz_all = 5;

% (1) Specify parameters by hand 
a0 = 0.5;
K = [0 35; 30 0];
Con = [18 16];
Coff = [1 1];
M_int = [0 1; -1 1];
hill = Inf;
noise = 0;
rcell = 0.2;
%cells = cells_hist{1};
lambda12 = 1.2;
lambda = [1 lambda12];
%p0 = s.p_ini;
p0 = [0.5 0.5];
%tmax =  s.tmax;
%gz = sqrt(N);
Rcell = rcell*a0;

% Settings
InitiateI = 0;
mcsteps = 0;
tmax = 10^4; % cut off simulation if t > tmax
%}
%% (2) Load parameters from saved trajectory
%
% with parameters saved as structure array 
% load data
data_folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\travelling waves';
%file = 'two_signal_mult_M_int0_1_-1_1_long_t_trans_wave_to_period10_trav_wave-v1';
%file = 'two_signal_mult_M_int0_1_-1_1_transient_wave_to_travelling_wave-v1';
%file = 'two_signal_mult_N225_initiateI1_I_ini_0p50_0p50_t_out_2436_period_15-v1';
file = 'two_signal_mult_M_int0_1_-1_1_hillInf_N225_a0_1p50_K0p0_22p0_11p0_4p0_Con18p0_16p0_Coff1p0_1p0-v1';
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
%cells = cells_hist{1};
lambda12 = s.lambda12;
lambda = [1 lambda12];
mcsteps = str2double(s.mcsteps);

p0 = s.p_ini;
%tmax =  s.tmax;
gz = sqrt(N);
Rcell = rcell*a0;
cell_type = zeros(N,1);

% Initial I
%
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
%}

% default file name
sim_ID = 'two_signal_mult';
I_ini_str = '';
if InitiateI
    I_ini_str = sprintf('_I_ini_%.2f_%.2f', I0(1), I0(2));
end
%}
%% temporary modification
K(1,2) = 9; 

%% First, calculate how many simulations are needed 
%%sim_to_do = zeros(numel(noise_all), 1);
% sim_to_do = zeros(numel(loopvar_all), 1);
sim_to_do = zeros(numel(hill_all), 1);

for idx1=1:numel(hill_all)
    %p0 = [p_all(idx1) p_all(idx2)];
    %noise = noise_all(idx1);
    hill = hill_all(idx1);

    % Count how many simulations have already been done
    %{
    %subfolder = strrep(sprintf('N%d strong int a0 %.1f', N, a0), '.', 'p');
    %folder = fullfile('L:\HY\Shared\Yiteng\two_signals\parameter set 2b', sprintf('N%d', N));
    %folder = fullfile('L:\HY\Shared\Yiteng\two_signals', 'sweep K22 new lattice', subfolder);
    %folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_Hill';
    %parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_noise\K12_9';
    %parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_Hill\K12_9';
    %parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_p0_set2';

    if exist(parent_folder, 'dir') ~= 7
        warning('Folder does not exist! ');
        break
    end
    %
    % subfolder
    subfolder = strrep(sprintf('ini_p1_%.2f_p2_%.2f', p0(1), p0(2)), '.', 'p');
    folder = fullfile(parent_folder, subfolder);
    if exist(folder, 'dir') ~= 7
        mkdir(folder);
    end
    %}
    folder = save_folder;
    
    % Filename pattern
    %pattern = strrep(sprintf('%s_N%d_p1_%.2f_p2_%.2f_K12_%d_t_out_%s_period_%s-v%s',...
    %        sim_ID, N, K(1,2), '(\d+)', '(\d+|Inf)',...
    %        '(\d+)'), '.', 'p');
    pattern = strrep(sprintf('%s_N%d_noise_%.2f_K12_%d_t_out_%s_period_%s-v%s',...
            sim_ID, N, noise, K(1,2), '(\d+)', '(\d+|Inf)',...
            '(\d+)'), '.', 'p');
    %pattern = strrep(sprintf('%s_N%d_hill_%.2f_K12_%d_t_out_%s_period_%s-v%s',...
    %        sim_ID, N, hill, K(1,2), '(\d+)', '(\d+|Inf)',...
    %        '(\d+)'), '.', 'p');
        
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

    %fprintf('N=%d, noise = %.2f sim to do: %d \n', N, noise, max_trials-filecount);
    fprintf('N=%d, hill = %.2f sim to do: %d \n', N, hill, max_trials-filecount);
    %fprintf('N=%d, p0 = [%.1f %.1f], sim to do: %d \n', N, p0(1), p0(2), max_trials-filecount);

    sim_to_do(idx1) = max_trials-filecount;
end
            
%% Then, do the simulations
for trial=1:max_trials
    for idx1=1:numel(hill_all)
        %noise = noise_all(idx1);
        hill = hill_all(idx1);
        %p0 = [p_all(idx1) p_all(idx2)];

        %fprintf('trial %d, N %d, noise %.2f \n', trial, N, noise);
        fprintf('trial %d, N %d, hill %.2f \n', trial, N, hill);
        %fprintf('trial %d, N %d, p0 = [%.1f %.1f] \n', trial, N, p0(1), p0(2));

        % skip simulation if enough simulations have been done
        if trial > sim_to_do(idx1)
            continue;
        end

        % ----------- simulation ------------------------------------
        [cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks_run2(...
            N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
            distances, positions, sim_ID, mcsteps, InitiateI, p0, I0, tmax, save_folder);
        %--------------------------------------------------------------------------
        %}
    end
end