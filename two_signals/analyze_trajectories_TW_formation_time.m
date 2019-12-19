%% Analyze TW formation time of saved trajectories across a range of parameters
clear all
close all
set(0, 'defaulttextinterpreter', 'tex');

%% Parameters
N = 225;
tmax = 10000;

%network = 33;
%num_params = 8; %
%nruns = 500; %number of runs per parameter set

network_all = [15 19 33 34 36];
num_params_all = [39 31 8 14 10];
nruns_all = [500 500 500 500 500];

% folder for saving figures
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_reliability';

%for network_idx = 1:numel(network_all)

%network = network_all(network_idx);
%num_params = num_params_all(network_idx);
%nruns = nruns_all(network_idx);
%% Load data files
%{
%%% PART I %%%
% folders for loading data
load_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_reliability';

subfolder = sprintf('TW_formation_network_%d', network);

names_all = cell(num_params, nruns);
for ii = 1:num_params
    subfolder2 = sprintf('Parameter_set_%d', ii);
    folder = fullfile(load_path, subfolder, subfolder2);

    listing = dir(folder);
    num_files = numel(listing)-2; %first two entries are not useful
    count = 0;
    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            count = count + 1;
            names_all{ii, count} = name;
        end
    end
end
%% Load data
filecount = zeros(num_params, 1); 
t_out_all = zeros(num_params, nruns); % final times
period_all = zeros(num_params, nruns); % periodicity test
t_onset_all = zeros(num_params, nruns); 
trav_wave_all = zeros(num_params, nruns); 
trav_wave_all_2 = zeros(num_params, nruns); 
hom_end_state_all = zeros(num_params, nruns); % whether end state is homogeneous

% --- Manually set ---
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_TW_params_%s_t_out_%s_period_%s-v%s',...
    N, '(\d+)', '(\d+)', '(\d+|Inf)', '(\d+)'); % '.' = anything
% --------------------

for ii=1:num_params
    subfolder = sprintf('TW_formation_network_%d', network);
    subfolder2 = sprintf('Parameter_set_%d', ii);
    folder = fullfile(load_path, subfolder, subfolder2);
    
    for jj=1:num_files
        name = names_all{ii, jj};
        if isempty(regexp(name, pattern, 'once')) % only load files matching a certain pattern
            continue
        else
            load( fullfile( folder, strcat(name, '.mat')), 'cells_hist',...
                'save_consts_struct', 'distances', 'period', 't_out', 't_onset');
            disp(name);
            %[tokens, ~] = regexp(name, pattern, 'tokens', 'match');

            filecount(ii) = filecount(ii) + 1;
            idx2 = filecount(ii);
            if idx2>nruns
                % only do up to nruns simulations
                continue
            end

            t_out_all(ii, idx2) = t_out;
            period_all(ii, idx2) = period;
            t_onset_all(ii, idx2) = t_onset;

            % test for TW
            a0 = save_consts_struct.a0;
            %digits = 3;
            [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                period, t_out, distances);
            trav_wave_all(ii, idx2) = trav_wave;
            trav_wave_all_2(ii, idx2) = trav_wave_2;

            % check end state
            cells_final = cells_hist{end};
            hom_end_state = size(unique(cells_final, 'rows'), 1)==1;
            hom_end_state_all(ii, idx2) = hom_end_state;
            %
        end
    end
end

% Check filecounts
disp(filecount);
%}
%% Save the loaded data
%{
subfolder = sprintf('TW_formation_network_%d', network);
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);

%
save_data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_reliability';
save( fullfile(save_data_path, strcat(fname_str, '.mat')),...
    'filecount', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
    'num_params', 'nruns', 'trav_wave_all', 'trav_wave_all_2',...
    'hom_end_state_all');
%}
%end

%% Plot averaged values: TW formation time
% (1) Load data
TW_frac_all = cell(5, 1);
TW_formation_time_all = cell(5, 1);
for network_idx = 1:5
    network = network_all(network_idx);
    num_params = num_params_all(network_idx);
    nruns = nruns_all(network_idx);

    subfolder = sprintf('TW_formation_network_%d', network);
    fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);

    % noise
    save_data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_reliability';
    load( fullfile(save_data_path, strcat(fname_str, '.mat') ),...
        'filecount', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
        'num_params', 'trav_wave_all', 'trav_wave_all_2', 'hom_end_state_all');
    
    % Classification: (1) TW end states, (2) oscillatory, (3) static, (4) infinite dynamics
    idx1 = (trav_wave_all_2);
    idx2 = (t_out_all < tmax & period_all < Inf & ~trav_wave_all_2);
    idx3 = (t_out_all < tmax & period_all == Inf );
    idx4 = (t_out_all == tmax);
    
    frac_all = [sum(idx1, 2) sum(idx2, 2)...
        sum(idx3, 2) sum(idx4, 2)]/nruns;
    
    TW_frac_all{network_idx} = frac_all(:, 1);
    
    % TW formation time
    TW_formation_time_all{network_idx} = t_onset_all(t_onset_all<Inf);
end

% Calculate mean and SEM
TW_frac_all_mean = zeros(5,1);
TW_frac_all_sem = zeros(5,1);
for i=1:5
    TW_frac_all_mean(i) = mean( TW_frac_all{i} );
    TW_frac_all_sem(i) = std( TW_frac_all{i} )/sqrt(numel(TW_frac_all{i}));
end

% Calculate mean and SEM of TW formation times
TW_formation_time_all_flatten = [];
boxplot_indices = [];
TW_formation_time_mean = zeros(5,1);
TW_formation_time_std = zeros(5,1);
TW_formation_time_sem = zeros(5,1);

edges = 0:250:10000;
TW_histcounts = zeros(5, numel(edges)-1);
for i=1:5
    TW_formation_time_mean(i) = mean(TW_formation_time_all{i});
    TW_formation_time_sem(i) = std(TW_formation_time_all{i})/sqrt(numel(TW_formation_time_all{i})) ;
    TW_formation_time_std(i) = std(TW_formation_time_all{i});
    TW_formation_time_all_flatten = [TW_formation_time_all_flatten; TW_formation_time_all{i}];
    boxplot_indices = [boxplot_indices; i*ones(numel(TW_formation_time_all{i}), 1)];
    TW_histcounts(i,:) = histcounts(TW_formation_time_all{i}, edges, 'Normalization', 'probability');
end

%% Plot average TW formation time per network
h = figure;
box on
hold on
boxplot(TW_formation_time_all_flatten, boxplot_indices)
%bar(TW_formation_time_mean);
%errorbar(1:5, TW_formation_time_mean, TW_formation_time_std, 'LineWidth', 2);
set(gca, 'XTick', 1:5, 'XTickLabels', sprintfc('%d', network_all));
xtickangle(45)
xlabel('Network');
ylabel('TW formation time');
set(gca, 'FontSize', 32);

qsave = 0;
fname = fullfile(save_path_fig, 'TW_formation_time', strcat('TW_formation_time_histogram',...
    sprintf('_nruns_%d_digits_%d', nruns, digits), '_all_networks_boxplot') );
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot all histograms together
% ------------> CONTINUE WORKING HERE <--------------------------------
h = figure;
b=bar3(TW_histcounts', 1);
set(b,'FaceAlpha', 0.5)
xlabel('Network');
set(gca, 'XTick', 1:5, 'XTickLabels', sprintfc('%d', network_all));
ylabel('TW formation time');
set(gca, 'YTick', 0.5:4:numel(edges)+0.5, 'YTickLabels', sprintfc('%d',...
    edges(1:4:end) ));
ylim([0 numel(edges)]);
zlabel('Probability');
set(gca, 'FontSize', 20);
pbaspect([1 2 1]);
view(-40, 20);
set(h, 'Units', 'inches', 'Position', [1 1 10 8]);

qsave = 0;
fname = fullfile(save_path_fig, 'TW_formation_time', strcat('TW_formation_time',...
    sprintf('_nruns_%d_digits_%d', nruns, digits), '_all_networks_histograms_together') );
save_figure(h, 10, 8, fname, '.pdf', qsave);
%--------------------------------------------------------------------
%% Plot TW formation time histograms per network
for i=1:5
    network = network_all(i);
    hist_data = TW_formation_time_all{i};
    
    edges = 0:250:10000;
    h=figure;
    hold on
    histogram(hist_data, edges, 'Normalization', 'probability');
    %set(gca, 'YScale', 'log');
    title(sprintf('%d simulations', numel(hist_data) ));
    xlabel('TW formation time');
    ylabel('Probability');
    ylim([0 0.2]);
    set(gca, 'FontSize', 32);
    set(h, 'Units', 'inches', 'Position', [1 1 10 8]);
    
    % test distribution    
    fitted_dist = fitdist(hist_data, 'Exponential');
    x = (edges(1:end-1)+edges(2:end))/2;
    pdf_calc = cdf(fitted_dist, edges(2:end)) - cdf(fitted_dist, edges(1:end-1));
    p = plot(x, pdf_calc, 'r-', 'LineWidth', 2);
    legend(p, sprintf('\\tau = %.0f \\pm %.0f', fitted_dist.mu, sqrt(fitted_dist.ParameterCovariance)));
    
    qsave = 0;
    fname = fullfile(save_path_fig, 'TW_formation_time', strcat('TW_formation_time_histogram',...
        sprintf('_nruns_%d_digits_%d_network_%d', nruns, digits, network), '_exp_fit'));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end


%% 
i=5;
hist_data = TW_formation_time_all{i};
fitted_dist = fitdist(hist_data, 'Exponential');
% x = (edges(2:end)+edges(1:end-1))/2; % evaluate distribution at bin centers

test_cdf = [hist_data, cdf('Exponential', hist_data, fitted_dist.mu)];
%h = kstest(hist_data, 'CDF', test_cdf)
[test_h, p, ksstat, cv] = kstest(hist_data, 'CDF', test_cdf);
disp(test_h);

% CONTROL: Test with randomly generated numbers from fitted distribution
%{
rand_nums = random('Exponential', fitted_dist.mu, 1000, 1);
test_cdf = [rand_nums, cdf('Exponential', rand_nums, fitted_dist.mu)];
[test_h, p, ksstat, cv] = kstest(rand_nums, 'CDF', test_cdf);
%}

figure;
hold on
cdfplot(hist_data)
x = linspace(0, tmax);
plot(x, cdf('Exponential', x, fitted_dist.mu) )
legend({'Data', 'Fitted distribution'});
