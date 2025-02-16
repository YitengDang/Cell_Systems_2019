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
p_all = 0:0.1:1; % loop variable (to be specified below)

% other settings
% p0=[0.5 0.5];
InitiateI = 0;
tmax = 10^4;

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

%p0 = s.p_ini;
%tmax =  s.tmax;
gz = sqrt(N);
Rcell = rcell*a0;
cell_type = zeros(N,1);

% Initial I
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

% default file name
sim_ID = 'two_signal_mult';
I_ini_str = '';
if InitiateI
    I_ini_str = sprintf('_I_ini_%.2f_%.2f', I0(1), I0(2));
end
%}

%% Temporary change
K(1,2) = 9;

%% First, calculate how many simulations are needed 
%sim_to_do = zeros(numel(noise_all), 1);
% sim_to_do = zeros(numel(loopvar_all), 1);
sim_to_do = zeros(numel(p_all));
for idx1=1:numel(p_all)
    for idx2=1:numel(p_all)
        p0 = [p_all(idx1) p_all(idx2)];
        %noise = noise_all(inner_idx);
        %hill = loopvar_all(inner_idx);
        
        % Count how many simulations have already been done
        %subfolder = strrep(sprintf('N%d strong int a0 %.1f', N, a0), '.', 'p');
        %folder = fullfile('L:\HY\Shared\Yiteng\two_signals\parameter set 2b', sprintf('N%d', N));
        %folder = fullfile('L:\HY\Shared\Yiteng\two_signals', 'sweep K22 new lattice', subfolder);
        %folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_Hill';
        parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_p0_set2';
        
        if exist(parent_folder, 'dir') ~= 7
            warning('Folder does not exist! ');
            break
        end
        subfolder = strrep(sprintf('ini_p1_%.2f_p2_%.2f', p0(1), p0(2)), '.', 'p');
        folder = fullfile(parent_folder, subfolder);
        if exist(folder, 'dir') ~= 7
            mkdir(folder);
        end
        
        % Filename pattern
        %pattern = strrep(sprintf('%s_N%d_initiateI%d%s_noise_%.2f_t_out_%s_period_%s%s-v%s',...
        %        sim_ID, N, InitiateI, I_ini_str, noise,...
        %        '(\d+)', '(\d+|Inf)', '\w*', '(\d+)'), '.', 'p');
        pattern = strrep(sprintf('%s_N%d_p1_%.2f_p2_%.2f_K12_%d_t_out_%s_period_%s-v%s',...
                sim_ID, N, p0(1), p0(2), K(1,2), '(\d+)', '(\d+|Inf)',...
                '(\d+)'), '.', 'p');
        
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
        fprintf('N=%d, p0 = [%.1f %.1f], sim to do: %d \n', N, p0(1), p0(2), max_trials-filecount);

        sim_to_do(idx1, idx2) = max_trials-filecount;
    end
end

%% Then, do the simulations
for trial=1:max_trials
    for idx1=1:numel(p_all) %numel(noise_all)
        for idx2=1:numel(p_all)
            %noise = noise_all(inner_idx);
            %hill = p_all(idx1);
            p0 = [p_all(idx1) p_all(idx2)];

            %fprintf('trial %d, N %d, noise %.2f \n', trial, N, noise);
            fprintf('trial %d, N %d, p0 = [%.1f %.1f] \n', trial, N, p0(1), p0(2));

            % skip simulation if enough simulations have been done
            if trial > sim_to_do(idx1, idx2)
                continue;
            end

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

            % always check within first t_ac time steps
            t_ac = 10^2; 
            while changed && period==Inf && t<t_ac
                %disp(t);
                %pause(0.2);
                t = t+1;
                cells = cellsOut;
                cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
                [period, t_onset] = periodicity_test_short(cells_hist); 
                %update_cell_figure_continuum(app, pos, dist, a0, cells, app.Time, cell_type, disp_mol, 0);
                [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
                    Rcell, Con, Coff, K, lambda, hill, noise);
            end
            % check periodically after t_ac time steps, with period t_check
            t_check = 10^3; 
            while changed && period==Inf && t<tmax
                %disp(t);
                %pause(0.2);
                t = t+1;
                cells = cellsOut;
                cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
                if mod(t, t_check)==0
                    [period, t_onset] = periodicity_test_short(cells_hist); 
                end
                %update_cell_figure_continuum(app, pos, dist, a0, cells, app.Time, cell_type, disp_mol, 0);
                [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
                    Rcell, Con, Coff, K, lambda, hill, noise);
            end
            %pause(1);

            t_out = t; %default t_out
            trav_wave = 0; % default trav_wave
            trav_wave_2 = 0;
            % if periodicity found, refine check to find period
            if period<Inf && t>t_ac
                [period, t_onset] = periodicity_test_detailed(cells_hist, t_check, period);
                t_out = t_onset + period;

                % also check if the solution is a traveling wave
                [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0, period, t_out); 
                % check over [t_out - period, t_out]
            end

            tmax_string = '';
            if changed && t==tmax
                tmax_string = '_tmax_reached';
            end

            fprintf('Final: t_out = %d, period %d \n', t_out, period);

            %--------------------------------------------------------------
            %% Save result
            % save folder
            subfolder = strrep(sprintf('ini_p1_%.2f_p2_%.2f', p0(1), p0(2)), '.', 'p');
            folder = fullfile(parent_folder, subfolder);
        
            % filename
            fname_str = strrep(sprintf('%s_N%d_p1_%.2f_p2_%.2f_K12_%d_t_out_%d_period_%s',...
                sim_ID, N, p0(1), p0(2), K(1,2), t_out, num2str(period)), '.', 'p');
            %fname_str = strrep(sprintf('%s_N%d_initiateI%d%s_randpos_mcsteps%d_K12_%d_t_out_%d_period_%s',...
            %    sim_ID, N, InitiateI, I_ini_str, mcsteps, K12, t_out, num2str(period)), '.', 'p');
            ext = '.mat';
            label = '';

            % check if filename already exists
            i=1;
            fname = fullfile(folder, strcat(fname_str, '-v', num2str(i), label, ext));
            while exist(fname, 'file') == 2
                i=i+1;
                fname = fullfile(folder, strcat(fname_str, '-v', num2str(i), label, ext));
            end
            
            % variables to be saved
            save_vars = {N, a0, K, Con, Coff, M_int, hill, noise, p0, rcell,...
                lambda12, sim_ID, I_ini_str, mcsteps};
            save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill',...
                'noise', 'p_ini', 'rcell', 'lambda12', 'sim_ID', ...
                'I_ini_str', 'mcsteps'};

            if InitiateI
                save_vars{end+1} = I0;
                save_vars_lbl{end+1} = 'I0';
            end

            save_consts_struct = cell2struct(save_vars, save_vars_lbl, 2);
            positions = pos;
            distances = dist;

            %
            save(fname, 'save_consts_struct', 'cells_hist', 't_out',...
                'changed', 'period', 't_onset', 'positions', 'distances',...
                'trav_wave', 'trav_wave_2');
            fprintf('Saved simulation: %s ; \n', fname);
            %}
            %--------------------------------------------------------------------------
        %}
        end
    end
end