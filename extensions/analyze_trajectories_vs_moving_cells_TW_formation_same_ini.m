%% Analyze saved trajectories across a range of parameters
clear all
close all
set(0, 'defaulttextinterpreter', 'tex');

%% Parameters
N = 225;
tmax = 10000;
sigma_D_all = [0 10.^[-4:0.5:-2]];
network = 15;
tmax = 2000;

var_all = sigma_D_all; %mcsteps_all;
label = 'network_15';
pset_all = 20; % parameter set

num_params = 1;
nruns = 100; %number of runs per parameter set

% folder for saving figures
save_path_fig = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW\TW_formation_network_15';

% folders for loading data
load_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
subfolder = sprintf('TW_formation_network_%d', network);

%% Get TW frequencies from loaded data
filecount = zeros(numel(var_all), num_params); 

% TW formation
periodicity_all = zeros( numel(sigma_D_all), num_params, nruns ); % whether each trajectory has at least one time-window in which it is periodic
periodicity_times = cell( numel(sigma_D_all), num_params, nruns ); % times at which solution is (temporarily) periodic
periodicity_periods = cell( numel(sigma_D_all), num_params, nruns ); % periods of the (temporarily) periodic times
TW_test_all = zeros( numel(sigma_D_all), num_params, nruns ); % whether each trajectory (transiently) forms a TW
TW_times_all = cell( numel(sigma_D_all), num_params, nruns ); % times of (Transient) TWs

% Load negative control (no noise) - fix later
%
load_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_reliability';
load_file = 'analyzed_data_TW_formation_network_15_nruns_500_digits_5';
%load_file = 'analyzed_data_TW_formation_network_15_nruns_500_digits_5_with_K_Con_data';
temp = load(fullfile(load_path, load_file));
%}

temp.idx_periodic = temp.period_all(pset_all, 1:nruns) < Inf;
periodicity_all(1, :, :) = temp.idx_periodic ;
%periodicity_periods(1, :, :) = num2cell(temp.period_all(1:num_params, 1:nruns)); % check whether correctly implemented!
TW_test_all(1, :, :) = temp.trav_wave_all_2(pset_all, 1:nruns);

%% Load data files
%
%%% PART I %%%
load_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';

% store parameter sets
Con_all_sim = zeros(num_params, 2); 
K_all_sim = zeros(num_params, 2, 2);

for idx_p_set = 1:numel(pset_all) %1:num_params
    disp(idx_p_set);
    p_set = pset_all(idx_p_set);
    subsubfolder = sprintf('Parameter_set_%d', p_set);
    folder = fullfile(load_path, subfolder, subsubfolder);
    
    % --- Manually set ---
    % mcsteps
    pattern = sprintf(...
        'two_signal_mult_N%d_sim_%d_params_%s_sigma_D_%s_tmax_%d-v%s',...
        N, p_set, '(\d+)', '(\dp\d+)', tmax,'(\d+)'); % '.' = anything
    
    % get all parameters from first simulation
    listing = dir(folder);
    num_files = numel(listing)-2; %first two entries are not useful
    for i = 1:num_files
        filename = listing(i+2).name;
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat') && ~isempty(regexp(filename, pattern, 'once'))
            load( fullfile( folder, filename), 'save_consts_struct');
            Con_all_sim(idx_p_set, :) = save_consts_struct.Con;
            K_all_sim(idx_p_set, :, :) = save_consts_struct.K;
            break
        end
    end
    
    % get filenames of all simulations
    names = {};
    listing = dir(folder);
    num_files = numel(listing)-2; %first two entries are not useful
    count = 0;
    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            count = count + 1;
            names{count} = name;
        end
    end
    
    % -------------------- 
    % load all simulations    
    for i=1:numel(names)
        if isempty(regexp(names{i}, pattern, 'once')) % only load files matching a certain pattern
            continue
        else
            %disp(names{i});
            load( fullfile( folder, strcat(names{i}, '.mat')), 'cells_hist',...
                'save_consts_struct', 'positions_all',...
                 'distances', 't_out', 'periodicity_vs_t', 'trav_wave_2_vs_t');

            [tokens, ~] = regexp(names{i}, pattern, 'tokens', 'match');

            % --- Manually set ---
            sigma_D = save_consts_struct.sigma_D;
            idx = find(sigma_D == sigma_D_all, 1);

            %noise = save_consts_struct.noise;
            %idx = find(noise == mcsteps_all, 1);
            % --------------------
            %{
            if num_params>1
                idx2 = str2double(tokens{1}{1});
            elseif num_params==1
                idx2 = 1;
            end        
            %}
            idx2 = idx_p_set;
            %disp(idx2);
            
            %
            if ~isempty(idx)
                filecount(idx, idx2) = filecount(idx, idx2) + 1;
                idx3 = filecount(idx, idx2);
                if idx3>nruns
                    % only do up to nruns simulations
                    continue
                end
                %disp(idx3);
                disp(names{i});
                
                %==================================================================
                % calculate periodicity_vs_t
                a0 = save_consts_struct.a0;
                [periodicity_vs_t, t_onset_over_time, trav_wave_2_vs_t] = ...
                    check_periodicity(cells_hist, a0, distances, tmax);

                % save new trajectory, also correct filename
                save_file_str = strrep(sprintf('two_signal_mult_N%d_sim_%d_params_%d_sigma_D_%.4f_tmax_2000-v1_new',...
                    N, p_set, idx3, sigma_D), '.', 'p');

                save( fullfile( folder, strcat(save_file_str, '.mat')), 'cells_hist',...
                    'save_consts_struct', 't_out', 'positions_all',...
                    'periodicity_vs_t', 't_onset_over_time', 'trav_wave_2_vs_t');
                %==================================================================
            
                % (2) for random initial state
                % N.B. This approach can only identify TWs that persist for times >= gz
                % is there any periodicity in the system?
                idx_temp = find(periodicity_vs_t>0);
                periodicity_all(idx, idx2, idx3) = ~isempty( idx_temp );
                fprintf('Periodic? %d \n', ~isempty(idx_temp));
                % store times of periodic solutions
                if periodicity_all(idx, idx2, idx3) 
                    periodicity_times{idx, idx2, idx3} = idx_temp;
                    periodicity_periods{idx, idx2, idx3} = periodicity_vs_t(idx_temp);
                end

                % (Transient) TWs
                idx_temp = find(trav_wave_2_vs_t == 1);
                TW_test_all(idx, idx2, idx3) = ~isempty( idx_temp );
                TW_times_all{idx, idx2, idx3} = idx_temp;
                %
            end
            %
        end
    end
    %}
end

% Save the loaded data
%
fname_str = sprintf('analyzed_data_%s_vs_sigma_D_num_params_%d_pset_20_nruns_%d_tmax_%d',...
    subfolder, num_params, nruns, tmax);
data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW\TW_formation_network_15';
save( fullfile(data_path, strcat(fname_str, '.mat')), 'sigma_D_all',...
    'filecount', 'num_params', 'nruns', 'tmax', 'periodicity_all',...
    'periodicity_times', 'periodicity_periods', 'TW_test_all', 'TW_times_all');

%% Load the saved data
fname_str = sprintf('analyzed_data_%s_vs_sigma_D_num_params_%d_pset_20_nruns_%d_tmax_%d',...
    subfolder, num_params, nruns, tmax);
data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW\TW_formation_network_15';
load( fullfile(data_path, strcat(fname_str, '.mat')), 'sigma_D_all',...
    'filecount', 'num_params', 'nruns', 'tmax', 'periodicity_all',...
    'periodicity_times', 'periodicity_periods', 'TW_test_all', 'TW_times_all');

%% Fraction of (transient) TWs (random initial state simulations)
count_TW = sum(sum(TW_test_all, 3), 2);
%y_data = count_TW/(nruns*num_params);
y_data = count_TW/(nruns*num_params);

x_data = sigma_D_all;
%x_data = sigma_D_all;
x_data(1) = sigma_D_all(2)/10; %[sigma_D_all(1)/10 sigma_D_all(2:end)];

%x_data = [sigma_D_all(1)/10 sigma_D_all];
%y_data = [TW_frac_neg_control(1); y_data];

h = figure;
hold on
plot(x_data, y_data, 'k--', 'LineWidth', 1.5 );
scatter(x_data, y_data, 'k', 'filled' );
set(gca, 'XScale', 'log');
box on

% Ticks
xticks = 10.^(-5:0);
% Tick labels
xtick_labels = sprintfc('10^{%d}', -5:0);
xtick_labels{1} = '0';
set(gca, 'FontSize', 32, 'XTick', xticks, 'XTickLabels', xtick_labels);

%{
h = figure;
hold on
bar( y_data, 'stacked');
set(gca, 'XTick', 1:numel(sigma_D_all), 'XTickLabel', sprintfc('%.3f', sigma_D_all) );
set(gca, 'YTick', 0:0.2:1);
%}
% plot errorbars
%sigma = std( TW_test_all(:,:), [], 2 );
% Plot confidence intervals
%{
z = 1.96;
CI_all = z*sigma./sqrt(nruns*num_params);
errorbars = CI_all;
%}
% plot standard deviation
%{
errorbars = sigma;
%}
%errorbar( bar_data, errorbars, '.', 'LineWidth', 3 )
%}
xlabel('Cell motility \sigma_D');
ylabel('Fraction TW');
set(gca, 'FontSize', 32);
set(h, 'Units', 'Inches', 'Position', [1 1 12 8]);
%ylim([0 0.08]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat(subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_vs_sigma_D_bar_plot_size_10_8_final'));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end
%% Detailed classification
% (to do)

%% Functions
function [periodicity_vs_t, t_onset_over_time, trav_wave_2_vs_t] = ...
    check_periodicity(cells_hist, a0, distances, tmax)
    periodicity_vs_t = []; % store all found periodicities over time
    t_onset_over_time = []; % store onset times of periodic trajectories
    %t_check = 1; % start checking periodicity from time 1
    trav_wave_2_vs_t = []; % whether simulation was a TW (according to second test, only constant p) for one period  
    

    for t=1:tmax
        pause(0.001);
        disp(t);
        %cells = cellsOut;
        %cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        %positions_all{end+1} = positions;

        % Check for periodicity
        [period, t_onset] = periodicity_test_short( cells_hist(1:t) );
        if period<Inf
            %[period, t_onset] = periodicity_test_detailed(cells_hist, t_check, period);
            [period, t_onset] = periodicity_test_short_reversed(cells_hist(1:t+1));
            periodicity_vs_t(t) = period;
        else
            periodicity_vs_t(t) = 0;
        end
        %periodicity_over_time(end+1) = period;
        t_onset_over_time(t) = t_onset;

        % Check for travelling waves
        if period<Inf
            cells_hist_temp = cells_hist(t-period+1:t+1);
            %dist = ones(N); % not needed, random value
            [~, trav_wave_2] = travelling_wave_test(cells_hist_temp, a0,...
                period, numel(cells_hist_temp)-1 , distances);
            trav_wave_2_vs_t(end+1) = trav_wave_2;
            fprintf('Time = %d, TW? %d \n', t, trav_wave_2);
        else
            trav_wave_2_vs_t(end+1) = 0;
        end
        %}       
    end
end