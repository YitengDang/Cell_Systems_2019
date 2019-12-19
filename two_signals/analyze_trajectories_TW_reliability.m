%% Analyze saved trajectories across a range of parameters
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

%% Process raw data
%{
for network_idx = 1:numel(network_all)

    network = network_all(network_idx);
    num_params = num_params_all(network_idx);
    nruns = nruns_all(network_idx);
    %% Load raw data files
    %
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
    Con_all = zeros(num_params, 2);
    K_all = zeros(num_params, 2, 2);

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

                % Store parameter data
                if jj==1 % first time: load parameters
                    Con_all(ii, :) = save_consts_struct.Con;
                    K_all(ii, :, :) = save_consts_struct.K;
                end

                % Store simulation data
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
    %
    %% Save the processed data
    %
    subfolder = sprintf('TW_formation_network_%d', network);
    fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5_with_K_Con_data', subfolder, nruns);

    %{
    save_data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_reliability';
    save( fullfile(save_data_path, strcat(fname_str, '.mat')),...
        'filecount', 'K_all', 'Con_all', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
        'num_params', 'nruns', 'trav_wave_all', 'trav_wave_all_2',...
        'hom_end_state_all');
    %}
end
%}
%% (1) Load data
TW_frac_all = cell(5, 1);
TW_formation_time_all = cell(5, 1);
Con_all_networks = cell(5, 1);
K_all_networks = cell(5, 1);

for network_idx = 1:5
    network = network_all(network_idx);
    num_params = num_params_all(network_idx);
    nruns = nruns_all(network_idx);

    subfolder = sprintf('TW_formation_network_%d', network);
    fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5_with_K_Con_data', subfolder, nruns);

    % reliability
    save_data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_reliability';
    load( fullfile(save_data_path, strcat(fname_str, '.mat') ),...
        'filecount', 'K_all', 'Con_all', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
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
    
    % Store parameter sets
    Con_all_networks{network_idx} = Con_all;
    K_all_networks{network_idx} = K_all;
end

% Calculate mean and SEM of reliability
TW_frac_all_mean = zeros(5,1);
TW_frac_all_sem = zeros(5,1);
for i=1:5
    TW_frac_all_mean(i) = mean( TW_frac_all{i} );
    TW_frac_all_sem(i) = std( TW_frac_all{i} )/sqrt(numel(TW_frac_all{i}));
end
%}
%% Plot average and std of TW reliability as bar graphs with errorbars
h=figure;
hold on
box on
bar(TW_frac_all_mean);
errorbar(TW_frac_all_mean, TW_frac_all_sem, 'LineWidth', 2)
ylim([0 1]);
%xlim([0.5 5.5]);
set(gca, 'XTick', 1:5, 'XTickLabels', sprintfc('%d', network_all));
xtickangle(45)
xlabel('Network');
ylabel('Reliability');
set(gca, 'FontSize', 32);

qsave = 0;
fname = fullfile(save_path_fig, strcat('analyzed_data_TW_formation',...
    sprintf('_nruns_%d_digits_%d', nruns, digits), '_classification_all_networks_bar_error'));
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot results as violin plot
% convert data into format for violin plot
TW_frac_converted = cell2mat(TW_frac_all);
TW_frac_cats = [];
for idx=1:numel(network_all)
    network = network_all(idx);
    num_params = num_params_all(idx);
    cats_temp = mat2cell(repmat(num2str(network), num_params, 1), num_params);
    TW_frac_cats = [TW_frac_cats; cats_temp{:}];
end
TW_frac_cats = categorical(cellstr(TW_frac_cats));
%%
h = figure;
violins = violinplot(TW_frac_converted, TW_frac_cats);
xlabel('Cellular Dialogue');
ylabel('Probability of forming wave');
set(gca, 'FontSize', 32);
box on
ylim([0 1]);

% save plot
qsave = 1;
fname = fullfile(save_path_fig, strcat('analyzed_data_TW_formation',...
    sprintf('_nruns_%d_digits_%d', nruns, digits), '_classification_all_networks_violin'));
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot detailed reliability figures
for network_idx=1:5 %[1 2 4 5]
    network = network_all(network_idx);
    num_params = num_params_all(network_idx);
    nruns = nruns_all(network_idx);
    Con_all = Con_all_networks{network_idx};
    K_all = K_all_networks{network_idx};
    
    %% Load the saved data
    %
    subfolder = sprintf('TW_formation_network_%d', network);
    fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);

    % noise
    save_data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_reliability';
    load( fullfile(save_data_path, strcat(fname_str, '.mat')),...
        'filecount', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
        'num_params', 'trav_wave_all', 'trav_wave_all_2', 'hom_end_state_all');

    %% 4-way classification
    % (1) TW end states, (2) oscillatory, (3) static, (4) infinite dynamics
    %
    idx1 = (trav_wave_all_2);
    idx2 = (t_out_all < tmax & period_all < Inf & ~trav_wave_all_2);
    idx3 = (t_out_all < tmax & period_all == Inf );
    idx4 = (t_out_all == tmax);
    % all(all(all(idx1 + idx2 + idx3 + idx4)))
    %}
    % Douwe's classification (4) static homogeneous, (3) oscillatory homogeneous, (1) TW end states, (2) infinite dynamics
    %{
    idx1 = (trav_wave_all_2);
    idx3 = (t_out_all < tmax & period_all == Inf & hom_end_state_all);
    idx4 = (t_out_all < tmax & period_all < Inf & ~trav_wave_all_2 & hom_end_state_all);
    idx2 = ones(size(idx1)) - (idx1+idx4+idx3); %(t_out_all == tmax);
    %all(all(all(idx1 + idx2 + idx3 + idx4)))
    %}
    % calculate fractions
    %frac_all = [sum(sum(idx1, 3), 2)/(num_params*nruns) sum(sum(idx2, 3), 2)/(num_params*nruns)...
    %    sum(sum(idx2, 3), 2)/(num_params*nruns) sum(sum(idx4, 3), 2)/(num_params*nruns)];
    frac_all = [sum(idx1, 2) sum(idx2, 2)...
        sum(idx3, 2) sum(idx4, 2)]/nruns;
    
    %% Plot fractions according to 4-way classification: unsorted
    h = figure;
    bar( frac_all, 'stacked');
    xlabel('Parameter set');
    ylabel('Fraction of simulations');
    ylim([0 1]);
    set(gca, 'FontSize', 32);
    box on
    % classificiation new
    % legend({'Travelling wave', 'Oscillatory', 'Static', 'Infinite dynamics'});
    % classification Douwe
    %legend({'Static homogeneous', 'Homogeneous oscillations', 'Pure wave', 'Infinite dynamics'});
    %legend({'Pure travelling wave', 'Infinite dynamics', 'Homogeneous, oscillatory', 'Homogeneous, static'});

    qsave = 0;
    fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_classification_bar_unsorted'));
    save_figure(h, 12, 8, fname, '.pdf', qsave);
    
    %% Plot fractions according to 4-way classification: sorted
    h = figure;
    [~, sort_idx] = sort(frac_all(:, 1), 'descend');
    frac_all_sorted = frac_all(sort_idx, :);
    % x_data = 1:num_params;

    % x-axis log-scale 
    %bar( log10(x_data), frac_all, 'stacked');
    %set(gca, 'XTick', -3:0, 'XTickLabels', sprintfc('10^{%d}', -3:0) );

    % x-axis evenly spread
    bar( frac_all_sorted, 'stacked');
    %set(gca, 'XTick', 1:numel(x_data), 'XTickLabels', sprintfc('%.3f', x_data) );
    xlabel('Parameter set');
    ylabel('Fraction of simulations');
    ylim([0 1]);
    set(gca, 'FontSize', 32);
    box on
    
    % classificiation new
    % legend({'Travelling wave', 'Oscillatory', 'Static', 'Infinite dynamics'});
    % classification Douwe
    %legend({'Static homogeneous', 'Homogeneous oscillations', 'Pure wave', 'Infinite dynamics'});
    %legend({'Pure travelling wave', 'Infinite dynamics', 'Homogeneous, oscillatory', 'Homogeneous, static'});

    qsave = 0;
    fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_classification_bar_sorted'));
    save_figure(h, 12, 8, fname, '.pdf', qsave);
    
    %% Plot fraction of TWs only
    %{
    %trav_wave_all_mean = sum(sum(trav_wave_all, 3), 2)/(num_params*nruns);
    trav_wave_all_2_mean = sum(trav_wave_all_2, 2)/(nruns);
    [~, sort_idx] = sort(trav_wave_all_2_mean, 'ascend');
    trav_wave_all_2_mean = trav_wave_all_2_mean(sort_idx, :);

    h = figure;
    x_data = 1:num_params;
    plot( x_data, trav_wave_all_2_mean, 'bo-', 'LineWidth', 1.5 );
    ylim([0 1]);
    %set(gca, 'XScale', 'log');

    % labels (change per case)
    xlabel('Parameter set');
    ylabel('Fraction TW');
    set(gca, 'FontSize', 32);

    qsave = 0;
    if qsave
        fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
            sprintf('_nruns_%d_digits_%d', nruns, digits), '_frac_TW_all_mean'));
        save_figure(h, 10, 8, fname, '.pdf', qsave);
    end
    %}
    
    %% Obtain TW formation times
    TW_formation_time_mean_all = zeros(num_params, 1);
    TW_formation_time_std_all = zeros(num_params, 1);
    TW_formation_time_count_all = zeros(num_params, 1);
    for i=1:num_params
        i_sorted = sort_idx(i);
        idx = t_onset_all(i_sorted,:)<Inf;
        TW_formation_time_mean_all(i) = mean(t_onset_all(i_sorted, idx));
        TW_formation_time_std_all(i) = std(t_onset_all(i_sorted, idx));
        TW_formation_time_count_all(i) = sum(idx);
    end
    
    h = figure;
    errorbar(1:num_params, TW_formation_time_mean_all,...
        TW_formation_time_std_all./TW_formation_time_count_all)
    
    %% Investigate parameter sets leading to TW formation
    % (1) K-Con maps
    %{
    % highest reliability 
    idx = sort_idx(1);
    Con_high = Con_all(idx, :);
    K_high = squeeze(K_all(idx, :, :));
    
    % lowest reliability
    idx = sort_idx(end);
    Con_low = Con_all(idx, :);
    K_low = squeeze(K_all(idx, :, :));
    %}
    
    h = figure;
    sgtitle(sprintf('Network %d', network), 'FontSize', 32);
    sz = 64; % marker size
    
    p11 = subplot(2,2,1);
    scatter( K_all(:, 1, 1), Con_all(:, 1), sz, frac_all(:, 1), 'filled' ); % 1<-1
    xlabel('K^{(11)}');
    ylabel('C_{ON}^{(1)}');
    c11=colorbar;
    ylabel(c11, 'TW reliability');
    
    p12 = subplot(2,2,2);
    scatter( K_all(:, 1, 2), Con_all(:, 2), sz, frac_all(:, 1), 'filled' ); % 1<-2
    xlabel('K^{(12)}');
    ylabel('C_{ON}^{(2)}');
    c12=colorbar;
    ylabel(c12, 'TW reliability');
     
    p21 = subplot(2,2,3);
    scatter( K_all(:, 2, 1), Con_all(:, 1), sz, frac_all(:, 1), 'filled' ); % 2<-1
    xlabel('K^{(21)}');
    ylabel('C_{ON}^{(1)}');
    c21=colorbar;
    ylabel(c21, 'TW reliability');
     
    set([p11 p12 p21], 'FontSize', 20);
    set([p11 p12 p21], 'xlim', [0 10^3], 'ylim', [0 10^3], 'clim', [0 1]);
    
    if network_idx>=3
        p22 = subplot(2,2,4);
        scatter( K_all(:, 2, 2), Con_all(:, 2), sz, frac_all(:, 1), 'filled' ); % 2<-2
        xlabel('K^{((22))}');
        ylabel('C_{ON}^{(2)}');
        c22=colorbar;
        ylabel(c22, 'TW reliability');
        
        set(p22, 'FontSize', 20);
        set(p22, 'xlim', [0 10^3], 'ylim', [0 10^3], 'clim', [0 1]);
    end
        
    qsave = 0;
    save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_reliability';
    fname_str = sprintf('Psets_network_%d_nruns_%d_by_reliability_K_Con_maps', network, nruns);
    fname = fullfile(save_folder, fname_str);
    save_figure(h, 18, 12, fname, '.pdf', qsave);
    %% Spider charts
    % Plot data as spider plot
    
    if network_idx<3
        K_idx = 1:3;
    else
        K_idx = 1:4;
    end
    
    %P_data = log10([Con_wave_sim, K_wave_sim(:,K_idx)]); % -> more generally, filter on M_int
    P_data = [Con_all, K_all(:, K_idx)];
    P_labels = {'C_{ON}^{(1)}', 'C_{ON}^{(2)}', 'K^{(11)}',...
        'K^{(21)}', 'K^{(12)}', 'K^{(22)}'};
    
    % for certain networks: swap molecules 1 <-> 2
    labelswap = '';
    if ~isempty(find(network_idx==[19 36], 1))
    %swap = 1;
    %if swap
        K_temp = zeros(size(K_all));
        K_temp(:,1,1) = K_all(:,2,2);
        K_temp(:,1,2) = K_all(:,2,1);
        K_temp(:,2,1) = K_all(:,1,2);
        K_temp(:,2,2) = K_all(:,1,1);
        P_data = [Con_wave_sim(:,2), Con_wave_sim(:,1), K_temp(:, K_idx)];
        labelswap = '_swapped_mol_1_2';
    end
    
    axes_interval = 3;
    
    % Plot spider plot
    h = figure;
    hold on
    
    TW_frac = TW_frac_all{network_idx};
    
    %color_all = repmat(TW_frac, 1, 3);
    %color_all =  'b'; 
    
    % Convert values in [0, 1] to colors through colormap
    %cmap = colormap('winter');
    cmap = colormap('jet');
    n_colors = size(cmap, 1);
    c_idx = floor(TW_frac./(1/n_colors));
    color_all = cmap(c_idx, :);
    
    % sort data
    [~, sort_idx] = sort( TW_frac, 'ascend');
    P_data = P_data(sort_idx, :);
    color_all = color_all(sort_idx, :);
    
    spider_plot_linear(P_data, P_labels([1:2 K_idx+2]), axes_interval,...
        'Color', color_all,...
        'Marker', 'o',...
        'LineStyle', '--',...
        'LineWidth', 1,...
        'MarkerSize', 10);
    
    c=colorbar('southoutside');
    set(c, 'FontSize', 20);
    xlabel(c, 'TW reliability');
    %title(sprintf('Network %d', network));
    %{
    spider_plot_linear_no_fill(P_data(sort_idx(1), :), P_labels([1:2 K_idx+2]), axes_interval, ...
        'Marker', 'none',...
        'LineStyle', '-',...
        'Color', [0 0 1],...
        'LineWidth', 3,...
        'MarkerSize', 2);
    %}
    
    qsave = 0;
    save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_reliability';
    fname_str = sprintf('Psets_network_%d_nruns_%d_by_reliability_spider_chart_v2_jet', network, nruns);
    fname = fullfile(save_folder, fname_str);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
    %% 
    %close all
end
%% fraction with period 15
%{
period_15_frac = sum( sum( period_all==15, 3 ), 2)/(num_params*nruns);
period_15_mult_frac = sum( sum( mod(period_all, 15)==0, 3 ), 2)/(num_params*nruns);

h = figure;
hold on
x_data = var_all;

plot( x_data, period_15_frac, 'bo-', 'LineWidth', 1.5  );
plot( x_data, period_15_mult_frac, 'ro-', 'LineWidth', 1.5  );
ylim([0 1]);
set(gca, 'FontSize', 20);
legend({'T=15', 'mod(T,15)=0'}, 'Location', 'sw');
set(gca, 'XScale', 'log');

% labels (change per case)
%xlabel('MC steps');
xlabel('Noise strength \alpha/K^{(ij)}');
ylabel('Fraction period 15');

% Tick labels 
%!!! Set manually !!!
%sel_idx = [1 6 11 16];
sel_idx = [1 3 5 7];
% last label
x_data(end+1) = 10^0;
A = sprintfc('10^{%d}', log10(x_data(sel_idx)) );
%A{end+1} = '10^{0}';
%{
A = cell( numel(mcsteps_all), 1 );
A{1} = '0';
A(2:end) = sprintfc('%d', mcsteps_all(2:end) );
%}
set(gca, 'FontSize', 32, 'XTick', x_data(sel_idx), 'XTickLabels', A);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_frac_period_15'));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end
%}
%% Analyze t_out
%{
% (1) <t_out>
%t_out_mean = mean(t_out_all, 3);

% (2) # trajectories reaching tmax
t_max_reached = (t_out_all == tmax);
t_max_count = sum(sum(t_max_reached, 3), 2);

h2 = figure;
x_data = var_all; 

plot( x_data, t_max_count/(num_params*nruns), 'bo-', 'LineWidth', 1.5  );
set(gca, 'FontSize', 32);
ylim([0 1]);
set(gca, 'XScale', 'log');

%xlabel('MC steps');
xlabel('Noise strength \alpha/K^{(ij)}');
ylabel('Fraction reaching t_{max}');
title(sprintf('t_{max}=10^{%d}', log10(tmax)));

%imagesc(p1, p2, t_max_count/nruns);
%set(gca, 'YDir', 'Normal');
%c = colorbar;
%caxis([0 1]);
%colormap('parula');
%xlabel('$$p_1$$')
%ylabel('$$p_2$$')
%ylabel(c, 'fraction');

% Tick labels 
%!!! Set manually !!!
%sel_idx = [1 6 11 16];
sel_idx = [1 3 5 7];
% last label
x_data(end+1) = 10^0;
A = sprintfc('10^{%d}', log10(x_data(sel_idx)) );
%A{end+1} = '10^{0}';
set(gca, 'FontSize', 32, 'XTick', x_data(sel_idx), 'XTickLabels', A);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_frac_reaching_tmax'));
    save_figure(h2, 10, 8, fname, '.pdf', qsave);
end
%}
%% TW detection algorithm test
%{
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\randomized lattice\TW_formation_network_15';
fname_str = 'two_signal_mult_N225_ini_state_rand_params_29_mcsteps_0_t_out_2707_period_225-v1';

fname = fullfile(folder, fname_str);
load(fname, 'cells_hist', 'save_consts_struct', 'distances', 'period', 't_out');

a0 = save_consts_struct.a0;
digits = 3;
[trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
    period, t_out, distances, digits);
%}
%%
% (3) <t_out> over trajectories that do not reach t_max
%{
t_out_mean_2 = zeros( numel(mcsteps_all), 1 );
for i=1:numel(t_out_all)
    if ~t_max_reached(i)
        [i1,i2,i3] = ind2sub(size(t_out_all), i); 
        t_out_mean_2(i1, i2) = t_out_mean_2(i1, i2) + t_out_all(i1,i2,i3);
    end
end
t_out_mean = t_out_mean_2./(nruns - t_max_count);

h3 = figure(3);

set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, '$$\langle t_{out} |  t_{out} < t_{max} \rangle$$', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
caxis([0 tmax]);
%ylim(c, [0 tmax]);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_mean_t_out_subset'));
    save_figure(h3, 10, 8, fname, '.pdf');
end
%}

%% Analyze periods
%{
% number of chaotic trajectories
n_chaotic = sum((period_all(:)==Inf).*(t_out_all(:) == tmax));

% data
period_data = period_all(period_all~=0);
uniq = unique(period_data);
C = categorical(period_data, uniq);

n_periodic = sum(period_data(:)~=Inf);

% distribution of periods
h4 = figure(4);
histogram(C);
xlabel('Period');
ylabel('Count');
title(sprintf('%d trajectories, %d periodic, %d chaotic',...
    sum(sum(filecount)), n_periodic, n_chaotic));
set(gca, 'FontSize', 20);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_periods_hist_chaotic'));
    save_figure(h4, 10, 8, fname, '.pdf');
end

%% average period (count only found periods)
period_mean = zeros(numel(p1), numel(p2));
period_count = zeros(numel(p1), numel(p2));
for i=1:numel(t_out_all)
    if period_all(i)>0 && period_all(i)<Inf
        [i1,i2,i3] = ind2sub(size(period_all), i); 
        period_count(i1, i2) = period_count(i1, i2) + 1;
        period_mean(i1, i2) = period_mean(i1, i2) + period_all(i1,i2,i3);
    end
end
period_mean = period_mean./period_count;

h5 = figure(5);
imagesc(p1, p2, period_mean);
set(gca, 'YDir', 'Normal');
c = colorbar;
caxis([0 max(uniq(uniq<Inf))])
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, 'Mean period (time steps)', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_mean_period'));
    save_figure(h5, 10, 8, fname, '.pdf');
end
%}

%% Fix data issue (K, Con values of parameter sets were missing)
%{
for network_idx = 1:5
network = network_all(network_idx);
num_params = num_params_all(network_idx);
nruns = nruns_all(network_idx);

% Load raw data files
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

Con_all = zeros(num_params, 2);
K_all = zeros(num_params, 2, 2);

% --- Manually set ---
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_TW_params_%s_t_out_%s_period_%s-v%s',...
    N, '(\d+)', '(\d+)', '(\d+|Inf)', '(\d+)'); % '.' = anything

for ii=1:num_params
    subfolder = sprintf('TW_formation_network_%d', network);
    subfolder2 = sprintf('Parameter_set_%d', ii);
    folder = fullfile(load_path, subfolder, subfolder2);
    
    for jj=1
        name = names_all{ii, jj};
        if isempty(regexp(name, pattern, 'once')) % only load files matching a certain pattern
            continue
        else
            load( fullfile( folder, strcat(name, '.mat')), 'cells_hist',...
                'save_consts_struct', 'distances', 'period', 't_out', 't_onset');
            disp(name);

            Con_all(ii, :) = save_consts_struct.Con;
            K_all(ii, :, :) = save_consts_struct.K;
        end
    end
end

%% Load and save
subfolder = sprintf('TW_formation_network_%d', network);
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);

save_data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_reliability';
load( fullfile(save_data_path, strcat(fname_str, '.mat')),...
    'filecount', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
    'num_params', 'nruns', 'trav_wave_all', 'trav_wave_all_2',...
    'hom_end_state_all');

save( fullfile(save_data_path, strcat(fname_str, '_with_K_Con_data.mat')),...
    'filecount', 'K_all', 'Con_all', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
    'num_params', 'nruns', 'trav_wave_all', 'trav_wave_all_2',...
    'hom_end_state_all');
%%
%{
network = 15;
nruns = 500;
subfolder = sprintf('TW_formation_network_%d', network);
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);

save_data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_reliability';
load( fullfile(save_data_path, strcat(fname_str, '_new.mat')),...
    'filecount', 'K_all', 'Con_all', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
    'num_params', 'nruns', 'trav_wave_all', 'trav_wave_all_2',...
    'hom_end_state_all');
%}
end
%}