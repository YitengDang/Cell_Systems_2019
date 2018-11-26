%% Obtain statistics by simulation for a fixed network
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'latex');
%% Parameters and settings
% Note: increasing nsim at n_pset is always possible. However, increasing
% n_pset leads to data sets that do not form a perfect LHS sample
n_pset = 10000; % number of parameter sets to do
nsim = 10; % number of simulations per parameter set

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

% network to examine
network_selected = [15 16 19 20 32 33 34 36 43]; % 14];

%%
for loop_idx=1:numel(network_selected)
    network = network_selected; %(loop_idx);
    disp(network);

    % folders
    save_folder_fig = fullfile('H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2',...
        sprintf('Network_%d', network) );

    % filename 
    fname_root = sprintf('Network_%d_analyzed_', network);
    %% Load data
    % Load analyzed data 
    folder = 'H:\My Documents\Multicellular automaton\data\two_signals\batch_sim_all_topologies';
    fname_str = 'batch_sim_analyzed_data_batch2';
    load(fullfile(folder,fname_str), 'network_all', 't_out_all', 'period_all',...
        'M_int_all', 'Con_all', 'K_all', 'tmax');

    % Get correct numbering
    network_orig = network_all(network); % original numbering of network

    % Select data 
    M_int = M_int_all{network_orig};
    t_out_all_temp = squeeze(t_out_all(network_orig, :, :));
    period_all_temp = squeeze(period_all(network_orig, :, :));
    Con_all_temp = squeeze(Con_all(network_orig, :, :));
    K_all_temp = squeeze(K_all(network_orig, :, :));

    disp(M_int);
    % Clear variables for memory
    clear Con_all
    clear K_all
    clear t_out_all
    clear period_all
    %% Load raw data
    %{
    subfolder1 = sprintf('Network_%d', network_orig);

    K_all = zeros(n_pset, 2, 2);
    Con_all = zeros(n_pset, 2);

    t_out_all = zeros(n_pset, nsim);
    period_all = zeros(n_pset, nsim);

    for idx1=1:n_pset
        subfolder2=sprintf('Param_%d', idx1);

        % get parameters
        fname = fullfile(parent_folder, subfolder1, subfolder2, 'parameters.mat');
        load(fname);
        disp(fname);

        K_all(idx1,:,:) = thisK;
        Con_all(idx1, :) = thisCon;

        for idx2=1:nsim
            % load file
            fname_str = sprintf('all_topologies_simulate-v%d.mat', idx2);
            fname = fullfile(parent_folder, subfolder1, subfolder2, fname_str);
            if exist(fname, 'file')~=2
                fname_str = sprintf('all_topologies_simulate-v%d_tmax_reached.mat', idx2);
                fname = fullfile(parent_folder, subfolder1, subfolder2, fname_str);
            end
            if exist(fname, 'file')~=2
                break
            end
            load(fname);
            %disp(fname);

            % save variables
            t_out_all(idx1, idx2) = t_out;
            period_all(idx1, idx2) = period;
        end
    end

    % Get interaction network from processed numbering (convert 1-44 to 1-81 numbering)
    % generate list of all networks
    done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
    M_int_all = {};
    network_all = [];
    for idx_n=1:3^4 
        [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], idx_n);
        M = [0 1 -1];
        M_int = [M(i11) M(i12); M(i21) M(i22)];
        if done(i11,i12,i21,i22)
            continue
        elseif idx_n==1
            continue
        else
            gM = [i22 i21; i12 i11];
            network_all(end+1) = idx_n;
            M_int_all{end+1} = M_int;
            done(i11,i12,i21,i22) = 1;
            done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
        end
    end

    % get current network
    disp(M_int_all{network});
    disp(network_all(network));

    %}
    %% 
    %{
    idx1 = 1:n_pset;
    figure;
    histogram(t_out_all_temp(idx1, :), 20);
    xlabel('t_{out}'); 
    ylabel('Count');
    %}
    %%
    %{
    idx1 = 1:n_pset;
    C = categorical(period_all_temp(idx1,:));
    figure;
    histogram(C);
    xlabel('Period');
    ylabel('Count');
    %}
    %% find long-lived dynamics in phase space
    % probably PCA would be better for visualization
    %{
    t_out_class = zeros(n_pset);
    for i=1:n_pset
        if mean(t_out_all_temp(i,:)) < 10
            t_out_class(i) = 1; % short-lived dynamics
        elseif mean(t_out_all_temp(i,:)) < 100
            t_out_class(i) = 2; % intermediate dynamics
        else 
            t_out_class(i) = 3; % long-lived dynamics
        end
    end

    idx1 = (t_out_class==1);
    idx2 = (t_out_class==2);
    idx3 = (t_out_class==3);
    %{
    figure;
    hold on
    plot(Con_all_temp(idx1,1), Con_all_temp(idx1,2), 'bo');
    plot(Con_all_temp(idx2,1), Con_all_temp(idx2,2), 'ro');
    plot(Con_all_temp(idx3,1), Con_all_temp(idx3,2), 'go');
    %}
    %% Find oscillation frequency for each parameter set
    osc_freq_all = zeros(n_pset, 1);
    for i=1:n_pset
        osc_freq_all(i) = sum(period_all_temp(i, :) < Inf)/nsim;
    end

    h1 = figure;
    plot(1:n_pset, osc_freq_all, 'bx');
    xlabel('Parameter set');
    ylabel('Oscillation frequency');
    set(gca, 'FontSize', 20);
    %% Plot unique periods per parameter set
    h2 = figure;
    hold on
    up_all = {};
    for i=1:n_pset
        up = unique(period_all_temp(i,:));
        up_all{i} = up;
        plot( i*ones(numel(up), 1), up, 'x');
    end
    pl = plot([1 n_pset], [4 4], 'r--');
    legend(pl, '4', 'Location', 'nw');
    xlabel('Parameter set');
    ylabel('Found periods');
    set(gca, 'FontSize', 20);
    %}
    %% Plot histogram of periods
    h2b = figure;
    C = categorical(period_all_temp);
    C_uniq = unique(C);
    idx = find(C~='Inf');
    histogram(C(idx), C_uniq(1:end-1))
    
    %% Plot t_out values per parameter set
    %{
    h3 = figure;
    hold on
    for i=1:n_pset
        plot( i*ones(nsim, 1), t_out_all_temp(i,:), 'x');
    end
    ylim([0 tmax]);
    xlabel('Parameter set');
    ylabel('$t_{out}$');
    set(gca, 'FontSize', 20);
    %}
    %% Find fraction of simulations with "complex dynamics" 
    % defined as either period>4 or t_out = tmax
    %{
    frac_complex_all = zeros(n_pset, 1);
    for i=1:n_pset
        idx1 = find(period_all_temp(i,:)>4 & period_all_temp(i,:)<Inf);
        idx2 = find(t_out_all_temp(i,:)==Inf);
        idx12 = union(idx1, idx2);
        if ~isempty(idx1)
            %fprintf('Parameter set: %d \n', i);
            %disp(period_all_temp(i,:));
            frac_complex_all(i) = numel(idx12)/nsim;
        end
    end

    % Plot fraction of parameter sets giving complex dynamics
    n_complex_psets = numel(find(frac_complex_all>0));
   
    h4 = figure;
    plot(1:n_pset, frac_complex_all, 'bx');
    ylim([0 1]);
    xlabel('Parameter set');
    ylabel('Fraction complex');
    title(sprintf('Complex dynamics in %d parameter sets', n_complex_psets));
    set(gca, 'FontSize', 20);
    %}
    %% Plot spider graph of parameter sets giving complex dynamics
    %{
    idx_complex = find(frac_complex_all>0);
    K_idx = 1:3;
    P_complex = log10([Con_all_temp(idx_complex,:) K_all_temp(idx_complex,K_idx)]);
    P_labels = {'$C_{ON}^{(1)}$', '$C_{ON}^{(2)}$', '$K^{(11)}$',...
        '$K^{(12)}$', '$K^{(21)}$'};
    axes_interval = 2;
    %axes_precision = 1;
    %
    
    spider_plot(P_complex, P_labels, axes_interval,...
        'Marker', 'o',...
        'LineStyle', '-',...
        'Color', [1 0 0],...
        'LineWidth', 0.5,...
        'MarkerSize', 2);
    %}
    %{
    hold on
    spider_plot([1 2 1 1 1 1; 1 2 1 1 1 2], P_labels, 1, 1,...
        'Marker', 'o',...
        'LineStyle', '-',...
        'LineWidth', 2,...
        'MarkerSize', 5);
    %
    % Title properties
    title(sprintf('Network %d', network),...
        'Fontweight', 'bold',...
        'FontSize', 16);
    set(gcf, 'Units', 'Inches', 'Position', [2 2 10 8]);
    h5 = gcf;
    
    % Save figure
    qsave = 0;
    save_figure(h5, 10, 8, fullfile(save_folder_fig, ...
        strcat(fname_root, 'parameters_complex_dynamics_spider')), '.pdf', qsave);
    %}
    %% Save figures
    %{
    qsave = 0;
    save_figure(h1, 10, 8, fullfile(save_folder_fig, ...
        strcat(fname_root, 'osc_freq_vs_param')), '.pdf', qsave);
    save_figure(h2, 10, 8, fullfile(save_folder_fig, ...
        strcat(fname_root, 'periods_vs_param')), '.pdf', qsave);
    save_figure(h3, 10, 8, fullfile(save_folder_fig, ...
        strcat(fname_root, 't_out_vs_param')), '.pdf', qsave);
    save_figure(h4, 10, 8, fullfile(save_folder_fig, ...
        strcat(fname_root, 'frac_complex_vs_param')), '.pdf', qsave);
    %}
    %%
    %close all
end