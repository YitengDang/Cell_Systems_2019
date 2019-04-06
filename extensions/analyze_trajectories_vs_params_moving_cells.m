%% Analyze saved trajectories across a range of parameters
clear all
close all
set(0, 'defaulttextinterpreter', 'tex');

%% Parameters
gz = 15;
N = gz^2;
tmax = 100;
%sigma_D_all = [0.001 0.003 0.005 0.01 0.03 0.05 0.1];
sigma_D_all = [0.001 10^(-2.75) 0.003 0.005 0.01 10^(-1.75) 0.03 0.05 0.1]; % 10^(-0.75)];
%sigma_D_all = [0.001 0.01 0.1];

num_params = 100;
%num_params = 30;
nruns = 3; %number of runs per parameter set
network = 15;

% Load data files
path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
subfolder = sprintf('TW_propagation_network_%d', network);

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
%% Reanalyze data and save (for reprocessing old data)
%{
path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
subfolder = 'temp_19';
load_folder = fullfile(path, subfolder);

path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
subfolder = 'temp_19_2';
save_folder = fullfile(path, subfolder);

pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_TW_params_%s_sigma_D_%s_tmax_%d-v%s',...
    N, '(\d+)', '\d+p\d+', tmax, '(\d+)'); 
%{
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_rand_params_%s_sigma_D_%s_tmax_%d-v%s',...
    N, '(\d+)', '\d+p\d+', tmax, '(\d+)'); 
%}
%fname_str = 'two_signal_mult_N225_ini_state_TW_params_1_sigma_D_0p001_tmax_100-v1';

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
        
        %------------------------------------------------------------------
        save(fullfile(save_folder, names{i}) ,...
            'cells_hist',...
            'positions_all',...
            'periodicity_vs_t',...
            'save_consts_struct',...
            't_out');
        %------------------------------------------------------------------
    end
end
%}

%% Load simulation data and analyze
%{
path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
subfolder = 'TW_formation_network_15';
folder = fullfile(path, subfolder);
%}
%-------Choice-------------------------------------------------------------
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_TW_params_%s_sigma_D_%s_tmax_%d-v%s',...
    N, '(\d+)', '\d+p\d+', tmax, '(\d+)'); % '.' = anything
%}
%{
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_rand_params_%s_sigma_D_%s_tmax_%d-v%s',...
    N, '(\d+)', '\d+p\d+', tmax, '(\d+)'); % '.' = anything#
%}
%--------------------------------------------------------------------------
names_ordered_all = cell( numel(sigma_D_all), num_params, nruns );

% data to store
filecount = zeros( numel(sigma_D_all), num_params ); 
TW_breaking_time_all = zeros( numel(sigma_D_all), num_params, nruns );
periodicity_trajectories = zeros( numel(sigma_D_all), num_params, nruns ); % whether each trajectory has at least one time-window in which it is periodic
periodicity_times = cell( numel(sigma_D_all), num_params, nruns ); % times at which solution is (temporarily) periodic
periodicity_periods = cell( numel(sigma_D_all), num_params, nruns ); % periods of the (temporarily) periodic times

for i=1:numel(names)
    if isempty(regexp(names{i}, pattern, 'once')) % only load files matching a certain pattern
        continue
    else
        disp(names{i});
        load( fullfile( folder, strcat(names{i}, '.mat')), 'cells_hist',...
            'save_consts_struct', 't_out',...
            'periodicity_vs_t');
        
        [tokens, ~] = regexp(names{i}, pattern, 'tokens', 'match');
        
        sigma_D = save_consts_struct.sigma_D;
        idx = find(sigma_D == sigma_D_all, 1);
        idx2 = str2double(tokens{1}{1});
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
            %
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
            %{
            % is there any periodicity in the system?
            idx_temp = find(periodicity_vs_t>0);
            periodicity_trajectories(idx, idx2, idx3) = ~isempty( idx_temp );
            fprintf('Periodic? %d \n', ~isempty(idx_temp));
            % store times of periodic solutions
            if periodicity_trajectories(idx, idx2, idx3) 
                periodicity_times{idx, idx2, idx3} = idx_temp;
                periodicity_periods{idx, idx2, idx3} = periodicity_vs_t(idx_temp);
            end
            %}
%------------END CHOICE--------------------------------------------------------
        end
        %}
    end
end
% reprocess
TW_breaking_time_all = TW_breaking_time_all + gz-1; % scale times into actual simulation times
TW_breaking_time_all(TW_breaking_time_all==gz) = 0; % TW broken before time gz

% Save analyzed data
%
% Save the loaded data
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);
save_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
save( fullfile(save_path, strcat(fname_str, '.mat')), 'sigma_D_all',...
    'filecount', 'tmax', 'num_params', 'nruns', 'TW_breaking_time_all' );

% folder for saving figures
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_formation_fixed_params_network_15\new'; %'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_moving_cells';
%}
%% Load analyzed data
%
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);
save_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\moving_cells_TW';
load( fullfile(save_path, strcat(fname_str, '.mat')), 'sigma_D_all',...
    'filecount', 'tmax', 'num_params', 'nruns', 'TW_breaking_time_all');

% folder for saving figures
%save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_formation_fixed_params_network_15\new'; %'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_moving_cells';
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_moving_cells';
%}
%% Plot persistence of TWs (bar graphs)
% Find breaking times of TWs (ini TW simulations)
h = figure;
count_Inf = sum(sum( TW_breaking_time_all == Inf, 3), 2);
count_0 = sum(sum( TW_breaking_time_all == 0, 3), 2);
count_other = num_params*nruns - count_Inf - count_0;

%bar_data = [count_Inf count_other count_0];
bar_data = [count_Inf count_other count_0]/(nruns*num_params);

% plot evenly spaced
%{
bar(bar_data, 'stacked');
set(gca, 'XTick', 1:numel(sigma_D_all), 'XTickLabel', sprintfc('%.3f', sigma_D_all) );
set(gca, 'YTick', 0:0.2:1);
%}
% plot log scale
%
b = bar(log10(sigma_D_all), bar_data, 'stacked', 'FaceColor', 'flat');
set(gca, 'XTick', -3:0, 'XTickLabel', sprintfc('10^{%d}', -3:0) );
set(gca, 'YTick', 0:0.2:1);
%}
%ylabel('Count');
ylabel('Fraction of simulations');
xlabel('Cell motility \sigma_D');
set(gca, 'FontSize', 32);
set(h, 'Units', 'Inches', 'Position', [1 1 12 8]);
legend({'Sustained', 'Broken at t>T', 'Broken at t\leqT' },...
    'Location', 'ne');

% play around with colors
cmap = colormap('parula');
colors = [cmap(1,:); cmap(ceil(size(cmap, 1)/2),:); cmap(size(cmap, 1) ,:)];
%colors = [1 0 0; 0 1 0; 0 0 1];
set(gca, 'colororder', colors);
for k = 1:size(bar_data,2)
    b(k).CData = repmat(colors(k,:), size(bar_data, 1), 1);
end

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_TW_preserved_vs_sigma_D_norm_prob'));
    save_figure(h, 12, 8, fname, '.pdf', qsave);
end
%% Plot persistence of TWs (fraction of simulations)
% 
h = figure;
% Plot only final time (t=1000, saturating value)
frac_sustained = sum(sum( TW_breaking_time_all == Inf, 3), 2)/(nruns*num_params);

% Plot multiple times
t1 = 100;
t2 = 50;
t3 = 30;
t4 = 10;
frac_sustained1 = sum(sum( TW_breaking_time_all > t1, 3), 2)/(nruns*num_params);
frac_sustained2 = sum(sum( TW_breaking_time_all > t2, 3), 2)/(nruns*num_params);
frac_sustained3 = sum(sum( TW_breaking_time_all > t3, 3), 2)/(nruns*num_params);
frac_sustained4 = sum(sum( TW_breaking_time_all > t4, 3), 2)/(nruns*num_params);
x_data = log10(sigma_D_all);
%y_data = [frac_sustained1 frac_sustained2 frac_sustained3 frac_sustained4];
y_data = frac_sustained;
plot(x_data, y_data, 'o-', 'LineWidth', 1.5);
set(gca, 'XTick', -3:0, 'XTickLabel', sprintfc('10^{%d}', -3:0) );
set(gca, 'YTick', 0:0.2:1);
ylim([0 1]);
ylabel('Fraction of sustained TWs');
xlabel('Cell motility \sigma_D');
set(gca, 'FontSize', 32);
set(h, 'Units', 'Inches', 'Position', [1 1 12 8]);
%legend(sprintfc('t=%d', [100 t2 t3 t4]));

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_frac_TW_preserved_vs_sigma_D_t1000_size_10_8'));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
    
    % save plot data
    save(fname, 'x_data', 'y_data');
end
%% Find filenames of specific simulations
%[i1, i2] = find( squeeze(TW_breaking_time_all(4, :, :) < Inf), 1);
%[i1, i2] = find( squeeze(TW_breaking_time_all(4, :, :) < Inf) & squeeze(TW_breaking_time_all(4, :, :) > 0), 1);
[i1, i2] = find( squeeze(TW_breaking_time_all(5, :, :) < Inf) & squeeze(TW_breaking_time_all(5, :, :) > 0));

soln = 1;
disp(names_ordered_all{4, i1(soln), i2(soln)})
TW_breaking_time_all(4, i1(soln), i2(soln))

%% Plot periodic times altogether
p_idx = 3;
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

%% Plot explicit data for given run
idx1 = 8;
idx2 = 1;
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