clear variables
close all
clc
set(0,'defaulttextinterpreter','latex');
%% Parameters and settings
% Settings
% Note: increasing nsim at n_pset is always possible. However, increasing
% n_pset leads to data sets that do not form a perfect LHS sample
n_pset = 100; % number of parameter sets to do
nsim = 100; % number of simulations per parameter set
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

% save folder
save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2';

%% Load saved data
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\batch_sim_all_topologies';
fname_str = sprintf('batch_sim_analyzed_data_batch2.mat');
load(fullfile(folder, fname_str));
n_networks = numel(M_int_all_reduced);

% filtered data
t_out_all_net = t_out_all(network_all, :, :);
period_all_net = period_all(network_all, :, :);

%% Calculate network distances (Hamming distance)
%graph_sim = zeros(n_networks);
graph_dist = zeros(n_networks);
for i=1:n_networks
    for j=1:n_networks
        %graph_dist(i,j) = 1;
        n_perm = 2;
        %M_sim_v1 = zeros(2); % measures how similar the two networks are
        M_dist = zeros(2); % measures how similar the two networks are
        for idx1=1:n_perm
            for idx2=1:n_perm
                M_int_1 = M_perm(M_int_all_reduced{i}, idx1);
                M_int_2 = M_perm(M_int_all_reduced{j}, idx2);
                %M_sim_v1(idx1, idx2) = sum(M_int_1(:)==M_int_2(:));
                M_dist(idx1, idx2) = sum(abs(M_int_1(:) - M_int_2(:)));
                %disp(M_int_1);
                %disp(M_int_2);
                %disp(M_sim(idx1, idx2));
            end
        end
        %graph_sim(i,j) = max(M_sim_v1(:));
        graph_dist(i,j) = min(M_dist(:));
    end
end
%graph_edges = (graph_sim==3); % If three links are the same, the networks are neighbours
graph_edges = (graph_dist==1); % If three links are the same, the networks are neighbours

%% Plot graphs
% Plot t_out max
t_out_max = max(t_out_all_net(:,:), [], 2);
h1 = figure;
G = graph(graph_edges);
% Layouts: circle, force(3), layered, subspace(3)
h = plot(G, 'Layout', 'force', 'MarkerSize', 20, 'LineWidth', 1.5, 'NodeCData',...
    t_out_max, 'NodeLabel', []);
for i=1:n_networks
   text(h.XData(i)+0.15,h.YData(i)-0.1, h.ZData(i), num2str(i), 'fontsize', 16, 'Color', 'k');
end
colormap(viridis);
c=colorbar;
c.FontSize = 20;
caxis([0 tmax]);
set(gca, 'FontSize', 20, 'XTick', [], 'YTick', []);
title('max($$t_{out}$$)')

% Plot t_out max
t_out_mean = mean(t_out_all_net(:,:), 2);
h2 = figure;
G = graph(graph_edges);
% Layouts: circle, force(3), layered, subspace(3)
h = plot(G, 'Layout', 'force', 'MarkerSize', 20, 'LineWidth', 1.5, 'NodeCData',...
    t_out_mean, 'NodeLabel', []);
for i=1:n_networks
	text(h.XData(i)+0.15,h.YData(i)-0.1, h.ZData(i), num2str(i), 'fontsize', 16, 'Color', 'k');
end
colormap(viridis);
c=colorbar;
c.FontSize = 20;
%caxis([0 1]);
set(gca, 'FontSize', 20, 'XTick', [], 'YTick', []);
title('$$\langle t_{out} \rangle$$')

% Save figures
qsave = 1;
fname_str = 'graph_plot_t_out_max';
save_figure(h1, 15, 10, fullfile(save_folder, fname_str),'.pdf', qsave);

qsave = 1;
fname_str = 'graph_plot_t_out_mean';
save_figure(h2, 15, 10, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Oscillation prevalence
frac_osc = zeros(n_networks, 1); % Fraction of parameter sets which are capable of generating oscillations
frac_osc_2 = zeros(n_networks, 1); % Net fraction of simulations which generate oscillations
for i=1:n_networks
    period_temp = squeeze(period_all_net(i,:,:));
    [idx1, idx2] = find(period_temp~=Inf);
    frac_osc(i) = numel(unique(idx1))/n_pset;
    frac_osc_2(i) = numel(idx1)/(n_pset*nsim);
end

% Plot oscillation prevalence
h3=figure;
G = graph(graph_edges);
% Layouts: circle, force(3), layered, subspace(3)
h = plot(G, 'Layout', 'force', 'MarkerSize', 20, 'LineWidth', 1.5, 'NodeCData',...
    frac_osc, 'NodeLabel', []);
for i=1:n_networks
   text(h.XData(i)+0.15,h.YData(i)-0.1, h.ZData(i), num2str(i), 'fontsize', 16, 'Color', 'k');
end
colormap(viridis);
c=colorbar;
c.FontSize = 20;
caxis([0 1]);
set(gca, 'FontSize', 20, 'XTick', [], 'YTick', [], 'ZTick', []);
title('Oscillation prevalence')

% Save figure
qsave = 1;
fname_str = 'graph_plot_osc_prevalence';
save_figure(h3, 15, 10, fullfile(save_folder, fname_str),'.pdf', qsave);
%% Test specific pair
%{
i = 5; 
j = 18;

n_perm = 2;
M_sim = zeros(2); % measures how similar the two networks are
for idx1=1:n_perm
    for idx2=1:n_perm
        M_int_1 = M_perm(M_int_all_reduced{i}, idx1);
        M_int_2 = M_perm(M_int_all_reduced{j}, idx2);
        M_sim(idx1, idx2) = sum(M_int_1(:)==M_int_2(:));
        disp(M_int_1);
        disp(M_int_2);
        disp(M_sim(idx1, idx2));
    end
end
this_sim = max(M_sim(:));
%}
%% Functions
M = [1 2; 3 4];
M_perm(M, 1)

function M_out = M_perm(M, perm)
    if perm==1
        M_out = M;
    elseif perm==2
        M_out = [M(2,2) M(2,1); M(1,2) M(1,1)];
    else
        warning('Wrong permutation');
        M_out = Inf;
    end
end