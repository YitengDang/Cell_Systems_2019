% Load parameter values from LHS sample
close all
clear all
maxNumCompThreads(6);
%warning off
%% Simulation parameters
max_trials = 100;

% Load LHS sample
path = 'H:\My Documents\Multicellular automaton\data\two_signals\LHS_sample';
nvals = 30;
v = 1;
fname_str = sprintf('LHS_values_K12_K21_Con1_Con2_nvals%d-v%d',...
        nvals, v);
load(fullfile(path, fname_str), 'K12_all', 'K21_all', 'Con1_all', 'Con2_all');

%% (1) input parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters
%Con = [18 16];
Coff = [1 1];
M_int = [0 1; -1 0];
%K = [3 12; 13 20]; % K(i,j): sensitivity of type i to type j molecules
lambda12 = 1.2;
lambda = [1 lambda12]; % diffusion length (normalize first to 1)
hill = Inf;
noise = 0;

% initial conditions
p0 = [0.5 0.5];
iniON = round(p0*N);
I0 = [0 0];
dI = 0.01;
InitiateI = 0; % 0: no, 1: yes

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1);

% pos, dist
[dist, pos] = init_dist_hex(gz, gz);

I_ini_str = '';
if InitiateI
    I_ini_str = sprintf('_I1_I2_%.2f_%.2f', I0(1), I0(2));
end
%}

%{
% TO DO: vectorize
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN1 = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(1)).*(lambda(1)./r)) ); % calculate signaling strength
fN2 = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(2)).*(lambda(2)./r)) ); % calculate signaling strength

% nearest neighbour interaction strength
fprintf('activator fij(a0) = %.4f \n', sinh(Rcell)*sum(exp((Rcell-a0)./lambda(1)).*(lambda(1)./a0)))
fprintf('inhibitor fij(a0) = %.4f \n', sinh(Rcell)*sum(exp((Rcell-a0)./lambda(2)).*(lambda(2)./a0)))
%}
%% (2) Load parameters from saved trajectory
%{
% with parameters saved as structure array 
% load data
%data_folder = 'H:\My Documents\Multicellular automaton\app\git_repository\release_2_1_full\data\time_evolution\sample_trajectories\typeIV';
data_folder = 'H:\My Documents\Multicellular automaton\app\git_repository\raw_current\data\time_evolution\parameter set 2b';
%data_folder = 'D:\Multicellularity\app\Multicellularity-2.1\data\time_evolution';
%file = 'two_signal_mult_M_int1_1_-1_-1_chaotic_state_tmax5000-v1.mat';
%file = 'two_signal_mult_M_int0_1_-1_1_long_t_trans_wave_to_period10_trav_wave-v1';
%file = 'two_signal_mult_M_int0_1_-1_1_transient_wave_to_travelling_wave-v1';
file = 'two_signal_mult_N225_initiateI1_I_ini_0p50_0p50_t_out_2436_period_15-v1';
%[file, data_folder] = uigetfile(fullfile(data_folder, '\*.mat'), 'Load saved simulation');
load(fullfile(data_folder, file));

s = save_consts_struct;
%N = s.N;
a0 = s.a0;
%K = s.K;
%Con = s.Con;
Coff = s.Coff;
M_int = s.M_int;
hill = s.hill;
noise = s.noise;
rcell = s.rcell;
cells = cells_hist{1};
lambda12 = s.lambda12;
lambda = [1 lambda12];
%p0 = s.p_ini;
p0 = [0.5 0.5];
%tmax =  s.tmax;
%gz = sqrt(N);
Rcell = rcell*a0;

% simulation parameters
%tmax = 10^4;
% Initial I

InitiateI = 0;
%{
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
%}

%% Loop
for idx=1:nvals
    % current loop variables
    K = [0 K12_all(idx); K21_all(idx) 0];
    Con = [Con1_all(idx) Con2_all(idx)];
    
    % save folder
    save_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\LHS_sample 1';

    %% Check existing files
    % Count how many simulations have already been done
    
    if exist(save_folder, 'dir') ~= 7
        warning('Folder does not exist! ');
    end
    sim_ID = 'two_signal_mult';

    % default file name
    I_ini_str = '';
    if InitiateI
        I_ini_str = sprintf('_I_ini_%.2f_%.2f', I0(1), I0(2));
    end

    filecount = 0;
    pattern = strrep(sprintf(...
            '%s_N%d_initiateI%d%s_K12_%.2f_K21_%.2f_Con1_%.2f_Con2_%.2f_t_out_%s_period_%s',...
            sim_ID, N, InitiateI, I_ini_str, K(1,2), K(2,1), Con(1), Con(2),...
            '(\d+)', '(\d+|Inf)'), '.', 'p');
        
    listing = dir(save_folder);
    num_files = numel(listing)-2;
    names = {};
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

    fprintf('idx = %d, sim to do: %d \n', idx, max_trials-filecount);
    %% Simulate
    for trial=1:max_trials-filecount
        fprintf('trial %d \n', trial);
        cells_hist = {};

        % generate initial lattice
        iniON = round(p0*N);
        cells = zeros(N, 2);
        for i=1:numel(iniON)
            cells(randperm(N,iniON(i)), i) = 1;
            if InitiateI && hill==Inf
                %fprintf('Generating lattice with I%d(t=0)... \n', i);
                dI = 0.1;
                [cells_temp, test, I_ini] = generate_I_new(cells(:, i), I0(i), I0(i)+dI, dist, a0);
                cells(:,i) = cells_temp;
                %fprintf('Generated initial I%d: %.2f; in range (Y=1/N=0)? %d; \n', i, I_ini, test);
            end
        end

        % store initial config
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        %-------------dynamics-----------------------------------------
        t = 0;
        period = Inf; %default values
        t_onset = Inf; 
        [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
                Rcell, Con, Coff, K, lambda, hill, noise);

        while changed && period==Inf %&& t<tmax
            t = t+1;
            cells = cellsOut;
            cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
            [period, t_onset] = periodicity_test_short(cells_hist); 
            [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
                Rcell, Con, Coff, K, lambda, hill, noise);
        end
        t_out = t; % save final time
        fprintf('t_out = %d, period %d \n', t_out, period);
        %--------------------------------------------------------------
        % Save result
        fname_str = strrep(sprintf(...
            '%s_N%d_initiateI%d%s_K12_%.2f_K21_%.2f_Con1_%.2f_Con2_%.2f_t_out_%d_period_%s',...
            sim_ID, N, InitiateI, I_ini_str, K(1,2), K(2,1), Con(1), Con(2),...
            t_out, num2str(period)), '.', 'p');
        ext = '.mat';
        label = '';

        % check if filename already exists
        i=1;
        fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), label, ext));
        while exist(fname, 'file') == 2
            i=i+1;
            fname = fullfile(save_folder, strcat(fname_str, '-v', num2str(i), label, ext));
        end

        if InitiateI
            save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, I0, rcell,...
                lambda12, sim_ID, I_ini_str};
            save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise', 'p_ini', 'I_ini', 'rcell',...
                'lambda12', 'sim_ID', 'I_ini_str'};
        else
            save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, rcell,...
                lambda12, sim_ID, I_ini_str};
            save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise', 'p_ini', 'rcell',...
                'lambda12', 'sim_ID', 'I_ini_str'};
        end

        save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);

        save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
            'changed', 'period', 't_onset');
            fprintf('Saved simulation: %s ; \n', fname);
        %--------------------------------------------------------------------------
    end
    %}
end