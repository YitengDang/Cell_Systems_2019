%% Analyze saved trajectories across a range of parameters
clear all
close all
set(0, 'defaulttextinterpreter', 'tex');

%% Parameters
gz = 15;
N = gz^2;
tmax = 1000;
sigma_D_all = [0 0.001 10^(-2.75) 0.003 0.005 0.01 10^(-1.75) 0.03 0.1 10^(-0.75)]; % 0.01 0.03 0.1]; %10.^[-2.75 -1.75 -0.75];
%sigma_D_all = [0.001 0.01 0.1];

num_params = 1;
nruns = 100; %number of runs per parameter set

% folder for saving figures
%save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_moving_cells';
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_formation_fixed_params_network_15';

%% Load data files (original files)
%{
path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
%subfolder = 'TW_propagation_network_19';
subfolder = 'TW_formation_network_15_fixed_params';

folder = fullfile(path, subfolder);
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
%}
%% Reanalyze data and save (for reprocessing old data)
%{
path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
subfolder = 'TW_formation_network_15_fixed_params';
load_folder = fullfile(path, subfolder);

path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
subfolder = 'TW_formation_network_15_fixed_params_reanalyzed';
save_folder = fullfile(path, subfolder);

pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_rand_fixed_params_sigma_D_%s_tmax_%d-v%s',...
    N, '\d+p\d+', tmax, '(\d+)'); 
pattern2 = sprintf(...
    'two_signal_mult_N%d_ini_state_rand_fixed_params_sigma_D_%s_tmax_%d',...
    N, '\d+p\d+', tmax);

%{
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_TW_params_%s_sigma_D_%s_tmax_%d-v%s',...
    N, '(\d+)', '\d+p\d+', tmax, '(\d+)'); 
%}
%{
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_rand_params_%s_sigma_D_%s_tmax_%d-v%s',...
    N, '(\d+)', '\d+p\d+', tmax, '(\d+)'); 
%}

for i=1:numel(names)
    if isempty(regexp(names{i}, pattern, 'once')) % only load files matching a certain pattern
        continue
    else
        disp(names{i});
        load( fullfile( load_folder, strcat(names{i}, '.mat')), 'cells_hist',...
            'save_consts_struct', 'positions_all', 'distances', 't_out',...
            'periodicity_over_time');
        
        %positions = positions_all{1};
        a0 = save_consts_struct.a0;
        cells_ini = cells_hist{1};
        p_ini = mean(cells_ini, 1);
        I_ini = zeros(1, 2);
        I_ini(1) = moranI(cells_ini(:,1), a0*distances);
        I_ini(2) = moranI(cells_ini(:,2), a0*distances);
        save_consts_struct.p_ini = p_ini;
        save_consts_struct.I_ini = I_ini;
        save_consts_struct.mcsteps = 0;
        %------------------------------------------------------------------
        % Check for periodicity
        periodicity_vs_t = zeros(tmax, 1);
        for tt=2:numel(cells_hist)
            this_period = periodicity_over_time(tt-1);
            if this_period<Inf
                t_check_init = 0;
                %[period, ~] = periodicity_test_detailed(cells_hist(1:tt), t_check_init,...
                %    this_period);
                [period, t_onset] = periodicity_test_short_reversed(cells_hist(1:tt));
                periodicity_vs_t(tt-1) = period;
            end
        end
        
        % Check for TW
        TW_times = find(mod(periodicity_vs_t, gz)==0 & periodicity_vs_t>0);
        trav_wave_2_vs_t = zeros(tmax, 1);
        for ii=1:numel(TW_times)
            this_time = TW_times(ii);
            this_period = periodicity_vs_t(this_time);

            cells_hist_temp = cells_hist(this_time-this_period+1:this_time+1);
            
            dist = ones(N); % not needed, random value
            [~, trav_wave_2] = travelling_wave_test(cells_hist_temp, a0,...
                this_period, numel(cells_hist_temp)-1 , dist);
            trav_wave_2_vs_t( this_time ) = trav_wave_2;

            fprintf('Time = %d, TW? %d \n', this_time, trav_wave_2);
            %
        end
        %------------------------------------------------------------------
        fname = fullfile(save_folder, strcat(fname_str, '.mat'));
        tokens = regexp(names{i}, pattern, 'tokens');
        v = str2double(tokens{1}{1});
        while exist(fname, 'file')==2
            v = v + 1;
            old_name = regexp(names{i}, pattern2, 'match');
            fname_str = sprintf('%s-v%d', old_name{1}, v);
            fname = fullfile(save_folder, strcat(fname_str, '.mat'));
        end

        save(fname,...
            'cells_hist',...
            'positions_all',...
            'periodicity_vs_t',...
            'save_consts_struct',...
            't_out',...
            'trav_wave_2_vs_t');
        %------------------------------------------------------------------
    end
end
%}

%% Load reanalyzed simulation data and analyze
% Load data files (original files)
%
path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
subfolder = 'TW_formation_network_15_fixed_params';
folder = fullfile(path, subfolder);

listing = dir(folder);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
names = {};
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
    end
end
%
%}
%-------Choice-------------------------------------------------------------
%{
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_TW_params_%s_sigma_D_%s_tmax_%d-v%s',...
    N, '(\d+)', '\d+p\d+', tmax, '(\d+)'); % '.' = anything
%}
%{
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_rand_params_%s_sigma_D_%s_tmax_%d-v%s',...
    N, '(\d+)', '\d+p\d+', tmax, '(\d+)'); % '.' = anything
%}
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_rand_fixed_params_sigma_D_%s_tmax_%d-v%s',...
    N, '\d+p\d+', tmax, '(\d+)');
pattern2 = sprintf(...
    'two_signal_mult_N%d_ini_state_rand_fixed_params_sigma_D_0p000_t_out_%s_period_%s',...
    N, '(\d+)', '(\d+|Inf)');
%--------------------------------------------------------------------------
names_ordered_all = cell( numel(sigma_D_all), num_params, nruns );

% data to store
filecount = zeros( numel(sigma_D_all), num_params ); 
% TW propagation
TW_breaking_time_all = zeros( numel(sigma_D_all), num_params, nruns );
% TW formation
periodicity_all = zeros( numel(sigma_D_all), num_params, nruns ); % whether each trajectory has at least one time-window in which it is periodic
periodicity_times = cell( numel(sigma_D_all), num_params, nruns ); % times at which solution is (temporarily) periodic
periodicity_periods = cell( numel(sigma_D_all), num_params, nruns ); % periods of the (temporarily) periodic times
TW_test_all = zeros( numel(sigma_D_all), num_params, nruns ); % whether each trajectory (transiently) forms a TW
TW_times_all = cell( numel(sigma_D_all), num_params, nruns ); % times of (Transient) TWs

for i=1:numel(names)
    if isempty(regexp(names{i}, pattern, 'once')) % only load files matching a certain pattern
        if isempty(regexp(names{i}, pattern2, 'once'))
            continue
        else
            % -- special case sigma_D = 0 --
            disp('sigmaD = 0');
            
            disp(names{i});
            load( fullfile( folder, strcat(names{i}, '.mat')), 'cells_hist',...
                'save_consts_struct', 'distances', 't_out', 'period');

            [tokens, ~] = regexp(names{i}, pattern, 'tokens', 'match');
        
            idx = 1;
            if num_params>1
                idx2 = str2double(tokens{1}{1});
            elseif num_params==1
                idx2 = 1;
            end
            
            filecount(idx, idx2) = filecount(idx, idx2) + 1;
            idx3 = filecount(idx, idx2);
            if idx3 > nruns
                % only do up to nruns simulations
                continue
            end
            %disp(idx3);
            
            % register file name        
            names_ordered_all{idx, idx2, idx3} = names{i};
            
            % TW times
            a0 = save_consts_struct.a0;
            digits = 3;
            [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                period, t_out, distances, digits);
            TW_test_all(idx, idx2, idx3) = trav_wave_2;
            % ------------------------------
            
        end
        
        %continue
    else
        disp(names{i});
        load( fullfile( folder, strcat(names{i}, '.mat')), 'cells_hist',...
            'save_consts_struct', 't_out',...
            'periodicity_vs_t', 'trav_wave_2_vs_t');
        
        [tokens, ~] = regexp(names{i}, pattern, 'tokens', 'match');
        
        sigma_D = save_consts_struct.sigma_D;
        idx = find(sigma_D == sigma_D_all, 1);
        if num_params>1
            idx2 = str2double(tokens{1}{1});
        elseif num_params==1
            idx2 = 1;
        end
        %disp(idx2);
        
        if ~isempty(idx)
            
            filecount(idx, idx2) = filecount(idx, idx2) + 1;
            idx3 = filecount(idx, idx2);
            if idx3 > nruns
                % only do up to nruns simulations
                continue
            end
            %disp(idx3);
            
            % register file name        
            names_ordered_all{idx, idx2, idx3} = names{i};
            
%------------CHOICE--------------------------------------------------------
            % (1) For initial state = TW
            %{
            % find first time the TW is broken (for initial wave)
            t_TW_broken = find(periodicity_vs_t(gz:end)~=gz, 1);
            if isempty(t_TW_broken)
                t_TW_broken = Inf; % TW never broken
            %elseif t_TW_broken==1 % TW broken before reaching 1 period
            end
            TW_breaking_time_all(idx, idx2, idx3) = t_TW_broken;
            %}
            
            % (2) for random initial state
            % N.B. This approach can only identify TWs that persist for times >= gz
            %
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
%------------END CHOICE--------------------------------------------------------
        end
        %
    end
end
% reprocess
TW_breaking_time_all = TW_breaking_time_all + gz-1; % scale times into actual simulation times
TW_breaking_time_all(TW_breaking_time_all==gz) = 0; % TW broken before time gz
periodicity_all = squeeze(periodicity_all);
%sum(periodicity_trajectories, 2)
%}
% Save analyzed data
%
% Save the loaded data
save_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';

% TW propagation
%{
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);
save( fullfile(save_path, strcat(fname_str, '.mat')), 'sigma_D_all',...
    'filecount', 'tmax', 'num_params', 'nruns', 'TW_breaking_time_all' );
%}

% TW formation
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5_with_sigmaD_0', subfolder, nruns);
save( fullfile(save_path, strcat(fname_str, '.mat')), 'sigma_D_all',...
    'filecount', 'tmax', 'num_params', 'nruns', 'periodicity_all',...
    'periodicity_times', 'periodicity_periods', 'TW_test_all', 'TW_times_all');
%}
%% Load analyzed data
%
save_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
subfolder = 'TW_formation_network_15_fixed_params';
% TW propagation
%{
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);
load( fullfile(save_path, strcat(fname_str, '.mat')), 'sigma_D_all',...
    'filecount', 'tmax', 'num_params', 'nruns', 'TW_breaking_time_all' );
%}

% TW formation
%fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);
%fname_str = sprintf('analyzed_data_TW_formation_network_15_fixed_params_reanalyzed_nruns_%d_digits_5', nruns);
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5_with_sigmaD_0', subfolder, nruns);

%load( fullfile(save_path, strcat(fname_str, '.mat')), 'sigma_D_all',...
%    'filecount', 'tmax', 'num_params', 'nruns', 'periodicity_all',...
%    'periodicity_times', 'periodicity_periods', 'TW_test_all', 'TW_times_all');
load( fullfile(save_path, strcat(fname_str, '.mat')), 'sigma_D_all',...
    'filecount', 'tmax', 'num_params', 'nruns', 'periodicity_all',...
    'periodicity_times', 'periodicity_periods', 'TW_test_all');
%}

%% Load the negative control data (other simulation set)
%{
nruns = 100;
subfolder = 'TW_formation_network_15_fixed_parameter_set';
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);

save_data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_with_noise';
load( fullfile(save_data_path, strcat(fname_str, '.mat')), 'trav_wave_all_2');
trav_wave_all_2_mean = sum(sum(trav_wave_all_2, 3), 2)/nruns;
TW_frac_neg_control = trav_wave_all_2_mean(1); % TW fraction of negative control (with sigma_D = 0), from other simulation set
%}

%% Fraction of (transient) TWs (random initial state simulations)
count_TW = sum(sum(TW_test_all, 3), 2);
%y_data = count_TW/(nruns*num_params);
y_data = count_TW(1:end-1)/(nruns*num_params);

x_data = sigma_D_all(1:end-1);
%x_data = sigma_D_all;
x_data(1) = sigma_D_all(2)/10; %[sigma_D_all(1)/10 sigma_D_all(2:end)];

%x_data = [sigma_D_all(1)/10 sigma_D_all];
%y_data = [TW_frac_neg_control(1); y_data];

h = figure;
plot(x_data, y_data, 'bo-', 'LineWidth', 1.5 );
set(gca, 'XScale', 'log');

% Ticks
xticks = 10.^(-4:0);
% Tick labels
xtick_labels = sprintfc('10^{%d}', -4:0);
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
%ylabel('Count');
%ylabel('Probability');
%ylabel('Frequency');
ylabel('Fraction TW');
set(gca, 'FontSize', 32);
set(h, 'Units', 'Inches', 'Position', [1 1 12 8]);
ylim([0 1]);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_TW_formation_vs_sigma_D_bar_plot'));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end
%% Find breaking times of TWs (ini TW simulations)
%{
h = figure;
count_Inf = sum(sum( TW_breaking_time_all == Inf, 3), 2);
count_0 = sum(sum( TW_breaking_time_all == 0, 3), 2);
count_other = num_params*nruns - count_Inf - count_0;

%bar_data = [count_Inf count_other count_0];
bar_data = [count_Inf count_other count_0]/(nruns*num_params);
%bar( sigma_D_all, [count_Inf count_other count_0], 'stacked');
bar( bar_data, 'stacked');
set(gca, 'XTick', 1:numel(sigma_D_all), 'XTickLabel', sprintfc('%.3f', sigma_D_all) );
set(gca, 'YTick', 0:0.2:1);
%ylabel('Count');
ylabel('Probability');
xlabel('$\sigma_D$');
set(gca, 'FontSize', 20);
set(h, 'Units', 'Inches', 'Position', [1 1 12 8]);
legend({'TW preserved', 'Broken at t>T', 'Broken at t\leqT' },...
    'Location', 'no');

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_TW_preserved_vs_sigma_D_norm_prob_v1b'));
    save_figure(h, 12, 8, fname, '.pdf', qsave);
end
%}
%% Analyze specific simulations
%{
%%%
%% Find filenames of specific simulations
%[i1, i2] = find( squeeze(TW_breaking_time_all(4, :, :) < Inf), 1);
%[i1, i2] = find( squeeze(TW_breaking_time_all(4, :, :) < Inf) & squeeze(TW_breaking_time_all(4, :, :) > 0), 1);
%[i1, i2] = find( squeeze(TW_breaking_time_all(5, :, :) < Inf) & squeeze(TW_breaking_time_all(5, :, :) > 0));

soln = 1;
name = names_ordered_all{4, i1(soln), i2(soln)};
disp(name)
TW_breaking_time_all(4, i1(soln), i2(soln))

%% Load trajectory
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW\TW_formation_network_15_fixed_params_reanalyzed';
fname_str = 'two_signal_mult_N225_ini_state_rand_fixed_params_sigma_D_0p001_tmax_1000-v10';
fname = fullfile(load_folder, fname_str);
load(fname);

a0 = save_consts_struct.a0;
dist = ones(N); % not required
% Check for TW
TW_times = find(mod(periodicity_vs_t, gz)==0 & periodicity_vs_t>0);
trav_wave_2_vs_t = zeros(tmax, 1);
for ii=1:numel(TW_times)
    this_time = TW_times(ii);
    this_period = periodicity_vs_t(this_time);
    
    cells_hist_temp = cells_hist(this_time-this_period+1:this_time+1);
    
    [~, trav_wave_2] = travelling_wave_test(cells_hist_temp, a0,...
        this_period, numel(cells_hist_temp)-1 , dist);
    trav_wave_2_vs_t( this_time ) = trav_wave_2;
    
    fprintf('Time = %d, TW? %d \n', this_time, trav_wave_2);
    %
end
%% Plot periodic times altogether
p_idx = 1;
plot_data = zeros(num_params*nruns, tmax+1);
temp = zeros(num_params*nruns,1);
for ii=1:num_params
    for jj=1:nruns
        %idx = sub2ind([num_params nruns], ii, jj);
        idx = (ii-1)*nruns + jj;
        temp(idx) = 1;
        %disp(idx);
        
        plot_data( idx, periodicity_times{p_idx, ii, jj} ) = periodicity_periods{p_idx, ii, jj};
    end
end

h = figure;
p_im = imagesc(plot_data);
set(p_im, 'AlphaData', plot_data>0);
colorbar;
%% Plot explicit data for given run
p_idx = 1;
idx1 = 1;
idx2 = 8;
disp(names_ordered_all{p_idx, idx1, idx2})

periodicity_vs_t = zeros(tmax+1, 1);
periodicity_vs_t(periodicity_times{p_idx, idx1, idx2}, 1) = periodicity_periods{p_idx, idx1, idx2};

h = figure;
plot(0:tmax, periodicity_vs_t, 'x')

%% Find possible TW solutions (periodicity ~ mult of gz)
% find onset time(s) of TW solutions
% calculate durations of possible TW solutions

%plot_data = zeros(num_params*nruns, tmax+1);
TW_possible = zeros(num_params*nruns,1);
for ii=1:num_params
    for jj=1:nruns
        %idx = sub2ind([num_params nruns], ii, jj);
        idx = (ii-1)*nruns + jj;
        temp(idx) = ~isempty(find(periodicity_periods{p_idx, ii, jj}==gz, 1));
        if temp(idx)
            disp(names_ordered_all{p_idx, ii, jj})
        end
    end
end
        
%% Periodicity vs time plot
%
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW\TW_propagation_network_19';
pset = 1;
run = 1;
%fname_str = 'two_signal_mult_N225_ini_state_TW_params_9_sigma_D_0p010_tmax_100-v2';
fname_str = sprintf('two_signal_mult_N225_ini_state_TW_params_%d_sigma_D_0p001_tmax_100-v%d', pset, run);
load( fullfile(load_folder, fname_str) );
%fname_str = 'two_signal_mult_N225_ini_state_TW_params_1_sigma_D_0p001_tmax_100-v1';
        
load( fullfile( load_folder, strcat(fname_str, '.mat')), 'cells_hist',...
    'save_consts_struct', 'positions_all', 't_out',...
    'periodicity_over_time');

%--------------------------------------------------------------------------
periodicity_vs_t_rev = Inf*ones(tmax, 1);
periodicity_vs_t_fwd = Inf*ones(tmax, 1);
for tt=1:numel(cells_hist)
    %[period, ~] = periodicity_test_detailed(cells_hist(1:tt), t_check_init,...
    %    this_period);
    [period, ~] = periodicity_test_short_reversed( cells_hist(1:tt) );
    periodicity_vs_t_rev(tt) = period;

    [period] = periodicity_test_short_forward( cells_hist(tt:end) );
    periodicity_vs_t_fwd(tt) = period;
end

h = figure;
hold on
plot( 0:tmax, periodicity_vs_t_rev );
plot( 0:tmax, periodicity_vs_t_fwd );
xlim([0 tmax]);
ymax = 5 + max( periodicity_vs_t_rev(periodicity_vs_t_rev<Inf));
ylim([0 ymax]);
xlabel('time');
ylabel('period');
set(gca, 'FontSize', 20);
%}