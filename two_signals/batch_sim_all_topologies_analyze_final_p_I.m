%% Obtain statistics on all possible topologies by simulation
% Examine final p, I
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'tex');

%% Parameters and settings
% Settings
% Note: increasing nsim at n_pset is always possible. However, increasing
% n_pset leads to data sets that do not form a perfect LHS sample
%n_pset = 10000; % number of parameter sets to do
%nsim = 10; % number of simulations per parameter set
tmax = 10000;

% Fixed parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];
hill = Inf;
noise = 0;
Coff = [1 1];
InitiateI = 0;

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% folders
load_folder = 'H:\My Documents\Multicellular automaton\data\two_signals\batch_sim_all_topologies\analysis_initial_p_I';
save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2';
subfolder1 = 'final_p';
subfolder2 = 'final_I';

%% Load saved data (first set: periodicity & t_out data)
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\batch_sim_all_topologies';
fname_str = sprintf('batch_sim_analyzed_data_batch2.mat');
load(fullfile(folder, fname_str), 't_out_all', 'period_all', 'M_int_all_reduced', 'network_all');
n_networks = numel(M_int_all_reduced);
n_pset = size(period_all, 2); % number of parameter sets
nsim = size(period_all, 3); % number of simulations per parameter set

% filtered data
t_out_all_net = t_out_all(network_all, :, :);
period_all_net = period_all(network_all, :, :);

clear t_out_all
clear period_all
%% 
% divide into subsets
%{
network_set{1} = [15  16	19	20	32	33	34	36	43]; % complex dynamics
network_set{2} = [2	5	8	10	11	12	13	14	18	21	22	23	25	27,...
                29	30	31	35	37	38	39	40	41	42	44]; % simple oscillations
network_set{3} = [1	 3	4	6	7	9	17	24	26	28]; % no oscillations
%}

% Get a list of all the networks
fname_str = 'batch_sim_analyzed_data_batch2';
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\batch_sim_all_topologies';
load(fullfile(folder, fname_str), 'network_all');

% 
networks_sel = [15  16	19	20	32	33	34	36	43]; %1:44;
for network=networks_sel
    %network = 1;
    fname_str = sprintf('batch_sim_analyzed_data_batch2_p_I_final_network_%d_old%d',...
        network, network_all(network) );
    fname = fullfile(load_folder, fname_str);
    disp(fname);
    load(fname, 'p_final_network', 'I_final_network');
    
    this_periods = period_all_net(network,:);
    %% Analyze final p
    %
    % Scatter plot
    h = figure;
    X = reshape(p_final_network(:, :, 1), n_pset*nsim, 1);
    Y = reshape(p_final_network(:, :, 2), n_pset*nsim, 1);
    scatter(X, Y, 'x');
    xlim([0 1]);
    ylim([0 1]);
    xlabel('Final $p^{(1)}$');
    ylabel('Final $p^{(2)}$');
    set(gca, 'FontSize', 24);

    qsave = 0;
    fname_str = sprintf('Final_p_scatter_network_%d', network);
    fname = fullfile(save_folder, subfolder1, fname_str);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
    
    %% Density plot (heat map)
    h = figure;
    Xedges = 0:0.05:1;
    Yedges = 0:0.05:1;
    Xcenters = (Xedges(1:end-1)+Xedges(2:end))/2;
    Ycenters = (Yedges(1:end-1)+Yedges(2:end))/2;
    [counts2, ~, ~] = histcounts2(X, Y, Xedges, Yedges,...
        'Normalization', 'probability');
    imagesc( Xcenters, Ycenters, log(counts2') );
    set(gca, 'YDir', 'normal');
    set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
    %xlim([0 1]);
    %ylim([0 1]);
    xlabel('Final $p^{(1)}$');
    ylabel('Final $p^{(2)}$');
    set(gca, 'FontSize', 24);
    c = colorbar;
    caxis([-10 0]); 
    set(c, 'YTick', -10:2:0, 'YTickLabel', sprintfc('10^{%d}', -10:2:0) );
    ylabel(c, 'Probability');
    scale = '_log'; % '_linear';

    qsave = 0;
    fname_str = sprintf('Final_p_density_map_network_%d%s', network, scale);
    fname = fullfile(save_folder, subfolder1, fname_str);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
    %% Analyze final I

    % Check maximum values
    fprintf('Min I = %.2f \n', min(I_final_network(:) ) );
    fprintf('Max I = %.2f \n', max(I_final_network(:) ) );

    % Scatter plot
    h = figure;
    X = reshape(I_final_network(:, :, 1), n_pset*nsim, 1);
    Y = reshape(I_final_network(:, :, 2), n_pset*nsim, 1);
    scatter(X, Y, 'x');
    xlim([-0.2 0.8]);
    ylim([-0.2 0.8]);
    xlabel('Final $I^{(1)}$');
    ylabel('Final $I^{(2)}$');
    set(gca, 'FontSize', 24);

    qsave = 0;
    fname_str = sprintf('Final_I_scatter_network_%d', network);
    fname = fullfile(save_folder, subfolder2, fname_str);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
    %% Density plot (heat map)
    h = figure;
    Xedges = -0.2:0.05:0.8;
    Yedges = -0.2:0.05:0.8;
    Xcenters = (Xedges(1:end-1)+Xedges(2:end))/2;
    Ycenters = (Yedges(1:end-1)+Yedges(2:end))/2;
    [counts2, ~, ~] = histcounts2(X,Y,Xedges,Yedges,...
        'Normalization', 'probability');
    imagesc( Xcenters, Ycenters, log(counts2') );
    set(gca, 'YDir', 'normal');
    set(gca, 'XTick', -0.2:0.2:0.8, 'YTick', -0.2:0.2:0.8);
    %xlim([0 1]);
    %ylim([0 1]);
    xlabel('Final $I^{(1)}$');
    ylabel('Final $I^{(2)}$');
    set(gca, 'FontSize', 32);
    c = colorbar;
    caxis([-10 0]); 
    set(c, 'YTick', -10:2:0, 'YTickLabel', sprintfc('10^{%d}', -10:2:0) );
    ylabel(c, 'Probability');
    scale = '_log'; % '_linear';

    qsave = 0;
    fname_str = sprintf('Final_I_density_map_network_%d%s_v2', network, scale);
    fname = fullfile(save_folder, subfolder2, fname_str);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
    %}
    %% Plot subset of points
    I_min = 0.3;
    I_1_final = I_final_network(:,:,1);
    I_2_final = I_final_network(:,:,2);
    p_1_final = p_final_network(:,:,1);
    p_2_final = p_final_network(:,:,2);
    
    % Scatter plot
    % select only periodic points
    %idx_subset = this_periods<Inf & this_periods>4;
    % select only periodic points with mod(T, gz) = 0
    idx_subset = (mod(this_periods, gz) == 0);
    % select spatially ordered points
    idx_subset2 = (I_1_final(:)>I_min & I_2_final(:)>I_min)';
    
    if sum(idx_subset)~=0
        num_in_gate = sum(idx_subset & idx_subset2);
        % Plot final p
        h = figure;
        hold on
        X = p_1_final(idx_subset);
        Y = p_2_final(idx_subset);
        scatter(X, Y, 'x');
        xlim([0 1]);
        ylim([0 1]);
        xlabel('Final $p^{(1)}$');
        ylabel('Final $p^{(2)}$');
        title(sprintf('$n_{ord} = %d, n_{tot}=%d$', num_in_gate, numel(X)));
        set(gca, 'FontSize', 32);
        box on
        
        qsave = 0;
        fname_str = strrep(sprintf(...
            'Final_p_scatter_complex_periods_network_%d_I_min_%.1f_v2',...
            network, I_min), '.', 'p');
        fname = fullfile(save_folder, subfolder1, fname_str);
        save_figure(h, 10, 8, fname, '.pdf', qsave);  
        
        % Plot final I
        h = figure;
        hold on
        X = I_1_final(idx_subset);
        Y = I_2_final(idx_subset);
        scatter(X, Y, 'o');
        plot([0.3 1], [0.3 0.3], 'r--');
        plot([0.3 0.3], [0.3 1], 'r--');
        plot([0.3 1], [1 1], 'r--');
        plot([1 1], [0.3 1], 'r--');
        xlim([-0.2 0.8]);
        ylim([-0.2 0.8]);
        xlabel('Final $I^{(1)}$');
        ylabel('Final $I^{(2)}$');
        title(sprintf('$n_{ord} = %d, n_{tot}=%d$', num_in_gate, numel(X)));
        set(gca, 'FontSize', 32);
        box on
        
        qsave = 1;
        fname_str = strrep(sprintf(...
            'Final_I_scatter_period_mod_gz_network_%d_I_min_%.1f_v2',...
            network, I_min), '.', 'p');
        fname = fullfile(save_folder, subfolder2, fname_str);
        save_figure(h, 10, 8, fname, '.pdf', qsave);   
    end
    %%
    close all
end

%% Estimate fraction of spatially ordered final states
% Arbitrary cutoff for telling when a pattern is spatially ordered
I_min = 0.3;

% Get a list of all the networks
fname_str = 'batch_sim_analyzed_data_batch2';
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\batch_sim_all_topologies';
load(fullfile(folder, fname_str), 'network_all');

% 
networks_sel = 1:44;
num_classes = 4; % number of different classes of trajectories to consider
num_ordered = zeros( numel(networks_sel), 1 );
I_mean_all = zeros( numel(networks_sel), num_classes, 2 );
I_std_all = zeros( numel(networks_sel), num_classes, 2 );
for network=networks_sel
    %network = 1;
    fname_str = sprintf('batch_sim_analyzed_data_batch2_p_I_final_network_%d_old%d',...
        network, network_all(network) );
    fname = fullfile(load_folder, fname_str);
    disp(fname);
    load(fname, 'I_final_network');
    
    % number of ordered networks
    num_ordered(network) = sum(sum(I_final_network(:,:,1)>I_min & I_final_network(:,:,2)>I_min));
    
    % final spatial order
    idx_class = {};
    
    % ordering 3
    idx_class{1} = squeeze(period_all_net(network, :, :) == Inf);
    idx_class{2} = squeeze(period_all_net(network, :, :) < 5);
    idx_class{3} = squeeze(period_all_net(network, :, :) < Inf & ...
        period_all_net(network, :, :) >= 5 & mod(period_all_net(network, :, :),gz)~=0);
    idx_class{4} = squeeze(period_all_net(network, :, :) < Inf & ...
        mod(period_all_net(network, :, :),gz)==0);
    
    % ordering 4
    %{
    idx_class{1} = squeeze(period_all_net(network, :, :) == Inf);
    idx_class{2} = squeeze(period_all_net(network, :, :) < 5);
    idx_class{3} = squeeze(period_all_net(network, :, :) < Inf & ...
        period_all_net(network, :, :) >= 5 & period_all_net(network, :, :) < gz);
    idx_class{4} = squeeze(period_all_net(network, :, :) < Inf & ...
        period_all_net(network, :, :) >= gz);
    %idx_class{4} = squeeze( mod(period_all_net(network, :, :),gz)==0 );
    %}
    for ii=1:num_classes
        idx = idx_class{ii};
        I1_data = I_final_network(:,:,1); I1_data = I1_data(idx); % flatten data
        I2_data = I_final_network(:,:,2); I2_data = I2_data(idx); % flatten data
        
        I_mean_all(network, ii, 1) = mean(I1_data);
        I_mean_all(network, ii, 2) = mean(I2_data);
        I_std_all(network, ii, 1) = std(I1_data);
        I_std_all(network, ii, 2) = std(I2_data); 
    end
    
end

% Plot fraction of ordered states
%{
h = figure;
bar(networks_sel, num_ordered/(n_pset*nsim));
% Plot as graph: see other file
%}

%% Plot average final I1, I2
mol = 1; % molecule #

% sort data
I_mean_all_sorted = I_mean_all; 
I_std_all_sorted = I_std_all;
idx_nan = isnan(I_mean_all_sorted);
I_mean_all_sorted(idx_nan) = -1; % set NaN values to -1 for ordering
sort_idx_final = 1:n_networks;
for ii = [1 2 3 4] %sort data one by one
    [~, sort_idx] = sort( I_mean_all_sorted(:, ii, mol), 'ascend' );
    I_mean_all_sorted(:, :, mol) = I_mean_all_sorted(sort_idx, :, mol);
    sort_idx_final = sort_idx_final(sort_idx);
end
I_std_all_sorted(:,:,mol) = I_std_all(sort_idx_final, :, mol);
I_mean_all_sorted( idx_nan(sort_idx_final, :, :) ) = NaN; % set the non-existing fractions to NaN again

h=figure;
hold on
for i=1:num_classes
    errorbar(I_mean_all_sorted(:,i,mol), I_std_all_sorted(:,i,mol), 'LineWidth', 2);
end
%legend({'T=\infty', 'T=2,4', 'T=3, T>5'}, 'Location', 'nw');
%legend({'T=\infty', 'T<5', 'T>=5, mod(T,gz)~=0', 'mod(T,gz)=0'}, 'Location', 'nw');
legend({'T=\infty', 'T<5', '5<=T<gz', 'T>=gz'}, 'Location', 'nw');

xlabel('Network');
ylabel(sprintf('Final I^{(%d)}', mol), 'Interpreter', 'tex');
set(gca, 'XTick', 1:n_networks, 'XTickLabel', sort_idx_final );
%set(gca, 'YTick', 0:0.2:1);
set(gca, 'FontSize', 20);
set(h, 'Units', 'Inches', 'Position', [1 1 24 6]);
box on
ylim([-0.2 0.8]);

% Save plot
qsave = 1;
save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2\final_order_by_network';
fname_str = sprintf('final_I%d_mean_by_network_errorbar_classification3', mol);
fname = fullfile(save_folder, fname_str);
save_figure(h, 24, 6, fname, '.pdf', qsave);
%% Plot distribution of final I per network
networks_sel = 15; %1:44;
for network=networks_sel
    %network = 1;
    fname_str = sprintf('batch_sim_analyzed_data_batch2_p_I_final_network_%d_old%d',...
        network, network_all(network) );
    fname = fullfile(load_folder, fname_str);
    disp(fname);
    load(fname, 'I_final_network');
    
    % Group by periodicity
    idx_class = {};
    idx_class{1} = squeeze(period_all_net(network, :, :) == Inf);
    idx_class{2} = squeeze(period_all_net(network, :, :) < 5);
    idx_class{3} = squeeze(period_all_net(network, :, :) < Inf & period_all_net(network, :, :) >= 5);
    %{
    idx_class{1} = squeeze(period_all_net(network, :, :) == Inf);
    idx_class{2} = squeeze(period_all_net(network, :, :) == 2);
    idx_class{3} = squeeze(period_all_net(network, :, :) == 4);
    idx_class{4} = squeeze(mod(period_all_net(network, :, :), gz)==0);
    idx_class{5} = squeeze(period_all_net(network, :, :) <Inf & period_all_net(network, :, :) ~=2 & ...
        period_all_net(network, :, :) ~=4 & mod(period_all_net(network, :, :), gz)~=0);
    %}
    
    I_mean = zeros(numel(idx_class), 2);
    I_std = zeros(numel(idx_class), 2);
    
    for ii=1:numel(idx_class)
        idx = idx_class{ii};
        I1_data = I_final_network(:,:,1); I1_data = I1_data(idx); % flatten data
        I2_data = I_final_network(:,:,2); I2_data = I2_data(idx); % flatten data
        
        I_mean(ii, 1) = mean(I1_data);
        I_mean(ii, 2) = mean(I2_data);
        I_std(ii, 1) = std(I1_data);
        I_std(ii, 2) = std(I2_data);  
    end
end
%% 
h = figure;
hold on
for i=1:numel(idx_class)
    errorbar(I_mean(i,1), I_mean(i,2),...
        I_std(i,1), I_std(i,1), I_std(i,2), I_std(i,2));
end
xlim([-0.2 0.8]);
ylim([-0.2 0.8]);
xlabel('I^{(1)}');
ylabel('I^{(2)}');
legend({'T=\infty', 'T=2,3,4', 'T>=5'});
%legend({'T=\infty', 'T=2', 'T=4', 'mod(T,gz)=0' 'Other T>1'});

%% Plot final I values for each period
mol = 1; %molecule

% fixed network
network = 43;
fname_str = sprintf('batch_sim_analyzed_data_batch2_p_I_final_network_%d_old%d',...
    network, network_all(network) );
fname = fullfile(load_folder, fname_str);
disp(fname);
load(fname, 'I_final_network');

% unique_periods = unique(period_all_net(network, :, :));
x_data = period_all_net(network, :, :);
y_data = I_final_network(:,:,1);
x_data(x_data==Inf) = 0.5;
h = figure;
scatter(log10(x_data(:)), y_data(:));

%% Histogram of final I for different conditions
mol = 2; %molecule

% fixed network
network = 14;
fname_str = sprintf('batch_sim_analyzed_data_batch2_p_I_final_network_%d_old%d',...
    network, network_all(network) );
fname = fullfile(load_folder, fname_str);
disp(fname);
load(fname, 'I_final_network');

idx = (period_all_net(network, :, :)==Inf);
I_data = I_final_network(:,:,mol);

figure;
nbins = 30;
histogram( I_data(idx), nbins );