% Generates all possible simulations by iterating over all initial states
% of a system
% -> Abandoned, does not work for even small systems (e.g. N=16)
close all
clear all
%% Simulation parameters
% grid size
gz = 4;
N = gz^2;

% folder to save simulations in
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\simulate_all_states_small_N';
subfolder = sprintf('N%d_%dstates', N, 4^N);
save_folder = fullfile(parent_folder, subfolder);

%% System parameters
%
% (1) Specify parameters by hand 
a0 = 1.5;
M_int = [0 1; -1 1];
K = [0 10; 11 4];
Con = [18 16];
Coff = [1 1];
hill = Inf;
noise = 0;
rcell = 0.2;
lambda12 = 1.2;
lambda = [1 lambda12];
p0 = [Inf Inf]; %[0.5 0.5];
I0 = [0 0];
Rcell = rcell*a0;

% get pos, dist
mcsteps = 0;
[positions, distances] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Settings
tmax = Inf; %10^4; % cut off simulation if t > tmax
%}
%% (2) Load parameters from saved trajectory
%{
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
%

I_ini_str = '';
if InitiateI
    I_ini_str = sprintf('_I_ini_%.2f_%.2f', I0(1), I0(2));
end
%}

%% Check how many simulations need to be done



%% Then, do the simulations
%{
for trial=1:sim_count
    for idx_loop=1:numel(N_all)
        N = N_all(idx_loop);

        subfolder = sprintf('N%d', N);
        save_folder = fullfile(parent_folder, subfolder);
        
        
            fprintf('N %d, trial %d \n', N, trial);
            %fprintf('trial %d, N %d, I0 = [%.2f %.2f] \n', trial, N, I0(1), I0(2));
            %fprintf('trial %d, N %d, noise %.2f \n', trial, N, noise);
            %fprintf('trial %d, N %d, p0 = [%.1f %.1f] \n', trial, N, p0(1), p0(2));

            % skip simulation if enough simulations have been done
            if trial > sim_to_do(idx_loop)
                disp('continue');
                continue;
            end

            [positions, distances] = initial_cells_random_markov_periodic(...
                sqrt(N), mcsteps, rcell);

            fname_str_template = sprintf('two_signal_mult_N%d', N);
            display_fig = 0;
            % ----------- simulation ------------------------------------
            [cells_hist, period, t_onset] = time_evolution_save_func_efficient_checks(...
                N, a0, Rcell, lambda, hill, noise, M_int, K, Con, Coff,...
                distances, positions, sim_ID, mcsteps, InitiateI, p0, I0, cells_ini, tmax,...
                save_folder, fname_str_template, display_fig);
            %--------------------------------------------------------------------------
            %}
        %end
    %end
%end