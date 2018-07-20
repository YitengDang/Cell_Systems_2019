% v3: use the randomization algorithm to place the cells on a
% different lattice
% v4: inner loop over K12 to keep the number of simulations with given
% parameters more or less constant
close all
clear all
maxNumCompThreads(6);
%warning off
%% Simulation parameters
max_trials = 100;

% lattice parameters
%gz = 15;
%N = gz^2;
gz_all = [25];

% loop over K12
K22_all = 5:5:25;

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
tmax = 10^5; % cut off simulation if t > tmax
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
K = s.K;
Con = s.Con;
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

%% First, calculate how many simulations are needed 
sim_to_do = zeros(numel(gz_all), numel(K22_all));
for gz_idx=1:numel(gz_all)  
    gz = gz_all(gz_idx);
    N = gz^2;
    for K22_idx=1:numel(K22_all)
        K22 = K22_all(K22_idx);
        % Count how many simulations have already been done
        subfolder = strrep(sprintf('N%d strong int a0 %.1f', N, a0), '.', 'p');
        %folder = fullfile('L:\HY\Shared\Yiteng\two_signals\parameter set 2b', sprintf('N%d', N));
        folder = fullfile('L:\HY\Shared\Yiteng\two_signals', 'sweep K22 new lattice', subfolder);

        if exist(folder, 'dir') ~= 7
            warning('Folder does not exist! ');
        end
        sim_ID = 'two_signal_mult';

        % default file name
        I_ini_str = '';
        if InitiateI
            I_ini_str = sprintf('_I_ini_%.2f_%.2f', I0(1), I0(2));
        end

        filecount = 0;
        pattern = strrep(sprintf('%s_N%d_initiateI%d%s_randpos_mcsteps%d_K22_%d_t_out_%s_period_%s%s-v%s',...
                sim_ID, N, InitiateI, I_ini_str, mcsteps, K22,...
                '(\d+)', '(\d+|Inf)', '\w*', '(\d+)'), '.', 'p');

        listing = dir(folder);
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

        fprintf('N=%d, K12 = %d, sim to do: %d \n', N, K22, max_trials-filecount);
        
        sim_to_do(gz_idx, K22_idx) = max_trials-filecount;
    end
end

%% Then, do the simulations
for trial=1:max_trials
    
    for gz_idx=1:numel(gz_all)    
        
        gz = gz_all(gz_idx);
        N = gz^2;

        cell_type = zeros(N,1);

        for K22_idx=1:numel(K22_all)
            K22 = K22_all(K22_idx);
            fprintf('trial %d, N %d, K12 %d \n', trial, N, K22);
            
            % skip simulation if enough simulations have been done
            if trial>sim_to_do(gz_idx, K22_idx)
                continue;
            end
            
            K(2,2) = K22;
            % ----------- simulation ------------------------------------
            cells_hist = {};

            % generate initial lattice
            %[dist, pos] = init_dist_hex(gz, gz);
            nodisplay = 1;
            [pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell, nodisplay);

            % generate initial state
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

            while changed && period==Inf && t<tmax
                %pause(0.2);
                t = t+1;
                cells = cellsOut;
                cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
                [period, t_onset] = periodicity_test_short(cells_hist); 

                %update_cell_figure_continuum(app, pos, dist, a0, cells, app.Time, cell_type, disp_mol, 0);
                [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
                    Rcell, Con, Coff, K, lambda, hill, noise);
            end
            t_out = t; % save final time
            if changed && t==tmax
                tmax_string = '_tmax_reached';
            else
                tmax_string = '';
            end
            %fprintf('t_out = %d, period %d \n', t_out, period);

            %--------------------------------------------------------------
            %% Save result
            fname_str = strrep(sprintf('%s_N%d_initiateI%d%s_randpos_mcsteps%d_K22_%d_t_out_%d_period_%s%s_temp',...
                sim_ID, N, InitiateI, I_ini_str, mcsteps, K22, t_out, num2str(period), tmax_string), '.', 'p');
            ext = '.mat';
            label = '';
            
            % check if filename already exists
            i=1;
            fname = fullfile(folder, strcat(fname_str, '-v', num2str(i), label, ext));
            while exist(fname, 'file') == 2
                i=i+1;
                fname = fullfile(folder, strcat(fname_str, '-v', num2str(i), label, ext));
            end

            if InitiateI
                save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, I0, rcell,...
                    lambda12, sim_ID, I_ini_str, mcsteps};
                save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise', 'p_ini', 'I_ini', 'rcell',...
                    'lambda12', 'sim_ID', 'I_ini_str', 'mcsteps'};
            else
                save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, rcell,...
                    lambda12, sim_ID, I_ini_str, mcsteps};
                save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill', 'noise', 'p_ini', 'rcell',...
                    'lambda12', 'sim_ID', 'I_ini_str', 'mcsteps'};
            end

            save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);
            positions = pos;
            distances = dist;
            save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
                'changed', 'period', 't_onset', 'positions', 'distances');
                fprintf('Saved simulation: %s ; \n', fname);
            %--------------------------------------------------------------------------
        end
    %}
    end
end