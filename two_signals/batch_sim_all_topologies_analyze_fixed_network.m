%% Obtain statistics by simulation for a fixed network
clear variables
close all
clc 
set(0, 'defaulttextinterpreter', 'tex');
%% Parameters and settings
% Note: increasing nsim at n_pset is always possible. However, increasing
% n_pset leads to data sets that do not form a perfect LHS sample
n_pset = 10000; % number of parameter sets to do
nsim = 10; % number of simulations per parameter set

% Fixed parameters
gz = 3;
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
%network_selected = [15 16 19 20 32 33 34 36 43]; % 14];
network = 15;

% folders
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_fixed_topology_vs_N';
%%
%for loop_idx=1:numel(network_selected)
    %network = network_selected; %(loop_idx);
    %disp(network);

    %===== folders ========================================================
    %subfolder1 = sprintf('Network_%d', network_orig);
    subfolder1 = sprintf('Network_%d_N%d', network, N);
    
    %save_folder_fig = fullfile('H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2',...
    %    sprintf('Network_%d', network) );
    save_folder_fig = fullfile(parent_folder); 
    
    % filename 
    fname_root = sprintf('Network_%d_N%d_analyzed_', network, N);
    %======================================================================
    %% Load analyzed data
    % =========== v1 ======================================================
    %{
    % Load analyzed data 
    load_folder = 'H:\My Documents\Multicellular automaton\data\two_signals\batch_sim_all_topologies';
    fname_str = 'batch_sim_analyzed_data_batch2';
    load(fullfile(load_folder, fname_str), 'network_all', 't_out_all', 'period_all',...
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
    %}
    % =========== batch_sim_fixed_topology ================================
    %{
    save_folder = parent_folder;
    fname_str = sprintf('analyzed_data_%s', subfolder1);
    load( fullfile(save_folder, fname_str), 'TW_strict_all', 'TW_loose_all', ...
        't_out_all', 'period_all', 'p_final_all', 'I_final_all',...
        'save_consts_struct', 'K_all', 'Con_all', 'fnames_all');
    %}
    % =====================================================================
    
    %% Load raw data
    %
    fnames_all = cell(n_pset, nsim);
    K_all = zeros(n_pset, 2, 2);
    Con_all = zeros(n_pset, 2);

    t_out_all = zeros(n_pset, nsim);
    period_all = zeros(n_pset, nsim);
    TW_strict_all = zeros( n_pset, nsim );
    TW_loose_all = zeros( n_pset, nsim );
    p_final_all = zeros( n_pset, nsim, 2 );
    I_final_all = zeros( n_pset, nsim, 2 );
    
    for idx1=1:n_pset
        subfolder2=sprintf('Param_%d', idx1);

        % get parameters
        path = fullfile(parent_folder, subfolder1, subfolder2);
        fname = fullfile(path, 'parameters.mat');
        load(fname);
        disp(fname);

        K_all(idx1,:,:) = thisK;
        Con_all(idx1, :) = thisCon;
        
        % ============ get filenames ======================================
        
        pattern = sprintf('Simulate_network_%d_params_%d_t_out_%s_period_%s-v%s',...
            network, idx1, '(\d+)', '(\d+|Inf)', '(\d+)');
        names = {};
        listing = dir( path );
        num_files = numel(listing)-2; %first two entries are not useful
        count = 0;
        for i = 1:num_files
            filename = listing(i+2).name;
            % remove extension and do not include txt files
            [~,name,ext] = fileparts(filename);
            if strcmp(ext, '.mat') && ~isempty(regexp(name, pattern, 'once'))
                count = count + 1;
                names{count} = name;
            end
        end

        % =================================================================
        for idx2=1:nsim
            % ============== specify filename =============================
            fname = fullfile(path, names{idx2});
            
            % load file
            %{
            fname_str = sprintf('all_topologies_simulate-v%d.mat', idx2);
            fname = fullfile(parent_folder, subfolder1, subfolder2, fname_str);
            if exist(fname, 'file')~=2
                fname_str = sprintf('all_topologies_simulate-v%d_tmax_reached.mat', idx2);
                fname = fullfile(parent_folder, subfolder1, subfolder2, fname_str);
            end
            if exist(fname, 'file')~=2
                break
            end
            %}
            %==============================================================
            
            load(fname, 't_out', 'period', 'save_consts_struct', 'cells_hist', 'distances');
            %disp(fname);
            
            % save variables
            fnames_all{idx1, idx2} = fname;
            t_out_all(idx1, idx2) = t_out;
            period_all(idx1, idx2) = period;
            
            cells_out = cells_hist{end};
            p_final_all(idx1, idx2, :) = sum(cells_out, 1)/N;
            I_final_all(idx1, idx2, 1) = moranI(cells_out(:, 1), a0*distances);
            I_final_all(idx1, idx2, 2) = moranI(cells_out(:, 2), a0*distances);
            
            % travelling wave test
            digits = -floor(log10(1/N)); % number of significant digits required to determine constancy of p
            a0 = save_consts_struct.a0;
            [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                period, numel(cells_hist)-1, distances, digits);
            
            TW_strict_all(idx1, idx2) = trav_wave;
            TW_loose_all(idx1, idx2) = trav_wave_2;
        end
    end
    %}
    %% Save analyzed data
    %
    save_folder = parent_folder;
    fname_str = sprintf('analyzed_data_%s', subfolder1);
    
    %{
    save_vars = {N, a0, K, Con, Coff, M_int, hill,...
        noise, p0, I0, rcell, Rcell, lambda12, lambda,...
        sim_ID, I_ini_str, tmax, mcsteps};
    save_vars_lbl = {'N', 'a0', 'K', 'Con', 'Coff', 'M_int', 'hill',...
        'noise', 'p0', 'I0', 'rcell', 'Rcell',  'lambda12', 'lambda', ...
        'sim_ID', 'I_ini_str', 'tmax', 'mcsteps'};
    parameters = cell2struct(save_vars, save_vars_lbl, 2);;
    %}
    
    save( fullfile(save_folder, fname_str), 'TW_strict_all', 'TW_loose_all', ...
        't_out_all', 'period_all', 'p_final_all', 'I_final_all',...
        'save_consts_struct', 'K_all', 'Con_all', 'fnames_all');
    %}
    %% Get interaction network from processed numbering (convert 1-44 to 1-81 numbering)
    %{
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
    %% plot t_eq (histogram)
    %
    %idx1 = 1:n_pset;
    %hist_data = t_out_all_temp(idx1, :);
    
    hist_data = t_out_all;
    
    C = categorical(hist_data);
    
    figure;
    histogram(C);
    xlabel('t_{out}'); 
    ylabel('Count');
    %}
    %% plot period (histogram)
    %
    %idx1 = 1:n_pset;
    %hist_data = period_all_temp(idx1,:);
    
    %hist_data = period_all;
    hist_data = period_all(period_all<Inf);
    
    C = categorical(hist_data);
    figure;
    histogram(C);
    xlabel('Period');
    ylabel('Count');
    %}
    
    %% Plot travelling waves
    % number of TWs
    num_TW_loose = sum(sum( TW_loose_all, 2 ), 1);
    num_TW_strict = sum(sum( TW_strict_all, 2 ), 1);
    
    % number of unique parameter sets 
    [idx1, ~] = find(TW_strict_all);
    numel( unique(idx1) )
    
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
    %{
    h2b = figure;
    C = categorical(period_all_temp);
    C_uniq = unique(C);
    idx = find(C~='Inf');
    histogram(C(idx), C_uniq(1:end-1))
    %}
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
    %% Plot spider graph of selected parameter sets
    %
    % select parameter sets
    %idx_select = find(frac_complex_all>0);
    [idx1, ~] = find(TW_strict_all);    
    %[idx1, ~] = find(TW_loose_all);  
    
    K_idx = 1:3;
    P_complex = [Con_all(idx1, :) K_all(idx1, K_idx)];
    %P_complex = log10([Con_all_temp(idx_select,:) K_all_temp(idx_select, K_idx)]);
    P_labels = {'C_{ON}^{(1)}', 'C_{ON}^{(2)}', 'K^{(11)}',...
        'K^{(12)}', 'K^{(21)}', 'K^{(22)}'};
    P_labels = P_labels( 1:2+numel(K_idx) );
    axes_interval = 4;
    range_values_in = repmat([0 1000], 5, 1);
    %
    
    h = figure;
    spider_plot(P_complex, P_labels, axes_interval, range_values_in,...
        'Marker', 'o',...
        'LineStyle', '-',...
        'Color', [1 0 0 0.3],...
        'LineWidth', 2,...
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
    %}
    
    % Save figure
    qsave = 1;
    fname_str = 'parameters_TW_strict_spider';
    %fname_str = 'parameters_TW_loose_spider';
    save_figure(h, 12, 12, fullfile(save_folder_fig, ...
        strcat(fname_root, fname_str)), '.pdf', qsave);
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
%end

%% Copy selected trajectories to folder
% load analyzed data
save_folder = parent_folder;
subfolder1 = sprintf('Network_%d_N%d', network, N);
fname_str = sprintf('analyzed_data_%s', subfolder1);
load( fullfile(save_folder, fname_str), 'TW_strict_all', 'TW_loose_all', ...
    't_out_all', 'period_all', 'p_final_all', 'I_final_all',...
    'save_consts_struct', 'K_all', 'Con_all', 'fnames_all');

% select trajectories
%[idx1, idx2] = find( TW_strict_all );
subfolder2 = 'batch_sim_selected_TW_strict';
idx = find( TW_strict_all );

%idx = find( TW_loose_all );
%subfolder2 = 'batch_sim_selected_TW_loose';

fname_sel = fnames_all(idx);

for ii=1:numel(fname_sel)
    fname = fname_sel(ii);
    fname = strcat(fname{1}, '.mat');
    
    [~, name, ext] = fileparts(fname);
    fname_new = fullfile(parent_folder, subfolder1, subfolder2,...
        strcat(name, ext));
    
    disp(fname_new);
    copyfile(fname, fname_new);
end