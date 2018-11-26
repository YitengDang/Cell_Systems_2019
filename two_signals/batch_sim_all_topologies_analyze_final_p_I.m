%% Obtain statistics on all possible topologies by simulation
% Examine final p, I
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'latex');

%% Parameters and settings
% Settings
% Note: increasing nsim at n_pset is always possible. However, increasing
% n_pset leads to data sets that do not form a perfect LHS sample
n_pset = 10000; % number of parameter sets to do
nsim = 10; % number of simulations per parameter set
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
networks_sel = 24:44;
for network=networks_sel
    %network = 1;
    fname_str = sprintf('batch_sim_analyzed_data_batch2_p_I_final_network_%d_old%d',...
        network, network_all(network) );
    fname = fullfile(load_folder, fname_str);
    disp(fname);
    load(fname);

    %% Analyze final p
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

    qsave = 1;
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

    qsave = 1;
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

    qsave = 1;
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
    set(gca, 'FontSize', 24);
    c = colorbar;
    caxis([-10 0]); 
    set(c, 'YTick', -10:2:0, 'YTickLabel', sprintfc('10^{%d}', -10:2:0) );
    ylabel(c, 'Probability');
    scale = '_log'; % '_linear';

    qsave = 1;
    fname_str = sprintf('Final_I_density_map_network_%d%s', network, scale);
    fname = fullfile(save_folder, subfolder2, fname_str);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
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
num_ordered = zeros( numel(networks_sel), 1 );
for network=networks_sel
    %network = 1;
    fname_str = sprintf('batch_sim_analyzed_data_batch2_p_I_final_network_%d_old%d',...
        network, network_all(network) );
    fname = fullfile(load_folder, fname_str);
    disp(fname);
    load(fname, 'I_final_network');
    
    num_ordered(network) = sum(sum(I_final_network(:,:,1)>I_min & I_final_network(:,:,2)>I_min));
end

% Plot fraction of ordered states

h = figure;
bar(networks_sel, num_ordered/(n_pset*nsim));

% Plot as graph: see other file
