%% Obtain statistics on all possible topologies by simulation, plot results as graphs
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'tex');
%% Parameters and settings
% Settings
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
save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2\graph_plots';

%% Load saved data
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

% divide into classes
class = cell(3, 1);
class{1} = [15	16	19	20  33	34	36	43]; % complex dynamics
class{2} = [2	5	8	10	11	12	13	14	18	21	22	23	25	27, ...
                29	30	31	32  35	37	38	39	40	41	42	44]; % simple oscillations
class{3} = [1	3	4	6	7	9	17	24	26	28]; % no oscillations

% Networks by number of repressive interactions
num_repressors = zeros(n_networks, 1);
num_repr_self = zeros(n_networks, 1); % self-interactions
for ii=1:n_networks
    num_repressors(ii) = sum(M_int_all_reduced{ii}(:)==-1);
    num_repr_self(ii) = sum(M_int_all_reduced{ii}([1 4])==-1);
end

%% Calculate network distances (Hamming distance)
%graph_sim = zeros(n_networks);
graph_dist = zeros(n_networks);
for ii=1:n_networks
    for j=1:n_networks
        %graph_dist(i,j) = 1;
        n_perm = 2;
        %M_sim_v1 = zeros(2); % measures how similar the two networks are
        M_dist = zeros(2); % measures how similar the two networks are
        for idx1=1:n_perm
            for idx2=1:n_perm
                M_int_1 = M_perm(M_int_all_reduced{ii}, idx1);
                M_int_2 = M_perm(M_int_all_reduced{j}, idx2);
                %M_sim_v1(idx1, idx2) = sum(M_int_1(:)==M_int_2(:));
                M_dist(idx1, idx2) = sum(abs(M_int_1(:) - M_int_2(:)));
                %disp(M_int_1);
                %disp(M_int_2);
                %disp(M_sim(idx1, idx2));
            end
        end
        %graph_sim(i,j) = max(M_sim_v1(:));
        graph_dist(ii,j) = min(M_dist(:));
    end
end
%graph_edges = (graph_sim==3); % If three links are the same, the networks are neighbours
graph_edges = (graph_dist==1); % If three links are the same, the networks are neighbours

%% Calculate x, y positions of networks
% get the number of interactions of each network
num_int = zeros(n_networks, 1); 
for ii=1:n_networks
    num_int(ii) = sum(sum( abs(M_int_all_reduced{ii}) ));
end

% number of networks with a given # of interactions
num_net_by_int = zeros(4, 1);
for ii=1:4
    num_net_by_int(ii) = sum(num_int == ii);
end

% find x and y positions of all networks, random arrangement (unused)
%{
x_graph = num_int;
y_graph = zeros(n_networks, 1);
count = zeros(4, 1);
for ii=1:n_networks
    this_num_int = num_int(ii); %number of interactions of this network
    y_offset = (num_net_by_int(this_num_int) - 1)/2; % readjust center

    y_graph(ii) = (count(this_num_int) - y_offset)/y_offset;
    count(this_num_int) = count(this_num_int) + 1;
end
%}

%% Calculate data for later usage
% tmax
t_out_max = max(t_out_all_net(:,:), [], 2);
t_out_max_log10 = log10(t_out_max);

% t_mean
t_out_mean = mean(t_out_all_net(:,:), 2);

% filter t_out on specific conditions
% (1) t_out avg over all with periodic ss
%{
[idx1, idx2] = find(period_all_net<Inf & period_all_net>4);
t_out_2 = zeros(size(t_out_all_net, 1), size(t_out_all_net, 2));
t_out_all_net_mat = t_out_all_net(:,:); % reshape into matrix
for ii=1:44
    idx_temp = find(idx1==ii);
    t_out_2(ii, idx2(idx_temp)) = t_out_all_net_mat(ii, idx2(idx_temp));
end
t_out_filtered = mean(t_out_2, 2);
%}

% (2) fraction with t_out > threshold
t_min = 50;
t_out_filtered = mean(t_out_all_net(:,:)>t_min, 2);

% Fraction oscillators
frac_osc = zeros(n_networks, 1); % Fraction of parameter sets which are capable of generating oscillations
frac_osc_2 = zeros(n_networks, 1); % Net fraction of simulations which generate oscillations
for ii=1:n_networks
    period_temp = squeeze(period_all_net(ii,:,:));
    [idx1, idx2] = find(period_temp~=Inf);
    frac_osc(ii) = numel(unique(idx1))/n_pset;
    frac_osc_2(ii) = numel(idx1)/(n_pset*nsim);
end
I_min = 0.3;

% Fraction spatially ordered
networks_sel = 1:44;
num_ordered = zeros( numel(networks_sel), 1 );
for network=networks_sel
    %network = 1;
    % Load files, find fraction of networks giving spatial order
    load_folder = 'H:\My Documents\Multicellular automaton\data\two_signals\batch_sim_all_topologies\analysis_initial_p_I';
    fname_str = sprintf('batch_sim_analyzed_data_batch2_p_I_final_network_%d_old%d',...
        network, network_all(network) );
    fname = fullfile(load_folder, fname_str);
    disp(fname);
    load(fname, 'I_final_network');
    
    num_ordered(network) = sum(sum(I_final_network(:,:,1)>I_min & I_final_network(:,:,2)>I_min));
end
frac_spatially_ordered = zeros(n_networks, 1);
frac_spatially_ordered(networks_sel) = num_ordered(networks_sel)/(n_pset*nsim);

%% Plot subset of networks as graph
subset = class{2}; 
graph_edges_temp = graph_edges(subset, subset);
networks_reduced = 1:n_networks;
networks_temp = networks_reduced(subset);

% get x, y positions
x_graph = num_int(subset);
y_graph = zeros(numel(x_graph), 1);
for ii=1:4
    idx = find(x_graph==ii);
    %disp(idx);
    idx_sorted = [];
    y_graph(idx) = linspace(-1, 1, numel(idx) );
end
%y_graph = y_graph(subset);
y_graph(1) = 0;

figure;
hold on
G = graph(graph_edges_temp);
% Layouts: circle, force(3), layered, subspace(3)
p = plot(G, 'Layout', 'force', 'MarkerSize', 0.1, 'LineWidth', 1.5, 'NodeLabel', []);
p.XData = x_graph;
p.YData = y_graph;
% plot text
for ii=1:numel(subset)
   text(p.XData(ii)+0.06,p.YData(ii)-0.05, p.ZData(ii), num2str(networks_temp(ii)),...
       'fontsize', 20, 'Color', 'k');
end

scatter(x_graph, y_graph, 20^2, t_out_max_log10(class{2}),...
    '^', 'filled');

%% (1) Plot max(t_out)
% find x and y positions, order by value of output variable
x_graph = num_int;
y_graph = zeros(n_networks, 1);
output_var = t_out_max;
for ii=1:4
    idx = find(num_int==ii);
    output_var_temp = output_var(idx);
    [output_var_sorted, sort_idx] = sort(output_var_temp, 'ascend');
    y_graph(idx(sort_idx)) = linspace(-1, 1, numel(idx) );
end

% Plot t_out max
h1 = figure;
hold on
G = graph(graph_edges);
% Layouts: circle, force(3), layered, subspace(3)
p = plot(G, 'Layout', 'force', 'MarkerSize', 0.1, 'LineWidth', 1.5, 'NodeCData',...
    log10(t_out_max), 'NodeLabel', []);

% Plot marker edges for classes
%scatter(x_graph(class{1}), y_graph(class{1}), p.MarkerSize^2, 'LineWidth', 2, 'MarkerEdgeColor', 'r');
%scatter(x_graph(class{2}), y_graph(class{2}), p.MarkerSize^2, 'LineWidth', 2, 'MarkerEdgeColor', 'g');
%scatter(x_graph(class{3}), y_graph(class{3}), p.MarkerSize^2, 'LineWidth', 2, 'MarkerEdgeColor', 'b');
scatter(x_graph(class{1}), y_graph(class{1}), 20^2, t_out_max_log10(class{1}),...
    'o', 'filled');
scatter(x_graph(class{2}), y_graph(class{2}), 20^2, t_out_max_log10(class{2}),...
    '^', 'filled');
scatter(x_graph(class{3}), y_graph(class{3}), 20^2, t_out_max_log10(class{3}),...
    's', 'filled');
p.XData = x_graph;
p.YData = y_graph;
for ii=1:n_networks
   text(p.XData(ii)+0.06,p.YData(ii)-0.05, p.ZData(ii), num2str(ii),...
       'fontsize', 20, 'Color', 'k', 'interpreter', 'tex');
end
%colormap(viridis);
colormap(jet);
c=colorbar;
c.FontSize = 24;
caxis([0 log10(tmax)]);
set(c, 'YTick', 0:log10(tmax), 'YTickLabel', sprintfc('10^{%d}',0:log10(tmax)) );
set(gca, 'FontSize', 24, 'XTick', 1:4, 'YTick', []);
xlabel('Number of interactions');
%ylabel(c, 'max($$t_{eq}$$)', 'Interpreter', 'tex')
title('Maximum equilibration time');
ylabel(c, 'Time steps', 'Interpreter', 'tex')

% Save figures
set(h1, 'Units', 'Inches', 'Position', [0.1 0.1 12 8]);
qsave = 0;
fname_str = 'graph_plot_t_out_max_log_rearranged_marked_v4_jet';
save_figure(h1, 12, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

%% (2) Plot t_out mean
% find x and y positions, order by value of output variable
x_graph = num_int;
y_graph = zeros(n_networks, 1);
output_var = t_out_mean;
for ii=1:4
    idx = find(num_int==ii);
    output_var_temp = output_var(idx);
    [output_var_sorted, sort_idx] = sort(output_var_temp, 'ascend');
    y_graph(idx(sort_idx)) = linspace(-1, 1, numel(idx) );
end

% Plot t_out mean
h2 = figure;
hold on
G = graph(graph_edges);
% Layouts: circle, force(3), layered, subspace(3)
p = plot(G, 'Layout', 'force', 'MarkerSize', 1, 'LineWidth', 1.5, 'NodeCData',...
    t_out_mean, 'NodeLabel', []);
% Plot marker edges for classes
scatter(x_graph(class{1}), y_graph(class{1}), 20^2, log10(output_var(class{1})),...
    'o', 'filled');
scatter(x_graph(class{2}), y_graph(class{2}), 20^2, log10(output_var(class{2})),...
    '^', 'filled');
scatter(x_graph(class{3}), y_graph(class{3}), 20^2, log10(output_var(class{3})),...
    's', 'filled');

p.XData = x_graph;
p.YData = y_graph;
for ii=1:n_networks
	text(p.XData(ii)+0.05,p.YData(ii)-0.05, p.ZData(ii), num2str(ii),...
        'fontsize', 16, 'Color', 'k', 'Interpreter', 'tex');
end
colormap(viridis);
%colormap(jet);
c=colorbar;
c.FontSize = 20;
%caxis([0 250]);
caxis([0 3]);
set(c, 'YTick', 0:3, 'YTickLabel', sprintfc('10^{%d}',0:3) );
set(gca, 'FontSize', 24, 'XTick', 1:4, 'YTick', []);
xlabel('Number of interactions');
%title('\langle t_{out} \rangle')
title('Mean equilibration time');
ylabel(c, 'Time steps', 'Interpreter', 'tex')

qsave = 0;
fname_str = 'graph_plot_t_out_mean_rearranged_marked_v4_viridis';
save_figure(h2, 12, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

%% (1,2 extra) Plot manually filtered t_out
% find x and y positions, order by value of output variable
x_graph = num_int;
y_graph = zeros(n_networks, 1);
output_var = t_out_filtered;
for ii=1:4
    idx = find(num_int==ii);
    output_var_temp = output_var(idx);
    [output_var_sorted, sort_idx] = sort(output_var_temp, 'ascend');
    y_graph(idx(sort_idx)) = linspace(-1, 1, numel(idx) );
end

% Plot t_out max
h = figure;
hold on
G = graph(graph_edges);
% Layouts: circle, force(3), layered, subspace(3)
p = plot(G, 'Layout', 'force', 'MarkerSize', 0.1, 'LineWidth', 1.5, 'NodeCData',...
    output_var, 'NodeLabel', []);

% Plot marker edges for classes
scatter(x_graph(class{1}), y_graph(class{1}), 20^2, output_var(class{1}),...
    'o', 'filled');
scatter(x_graph(class{2}), y_graph(class{2}), 20^2, output_var(class{2}),...
    '^', 'filled');
scatter(x_graph(class{3}), y_graph(class{3}), 20^2, output_var(class{3}),...
    's', 'filled');
p.XData = x_graph;
p.YData = y_graph;
for ii=1:n_networks
   text(p.XData(ii)+0.06,p.YData(ii)-0.05, p.ZData(ii), num2str(ii),...
       'fontsize', 20, 'Color', 'k', 'interpreter', 'tex');
end
colormap(viridis);
%colormap(jet);
c=colorbar;
c.FontSize = 24;
%caxis([0 log10(tmax)]);
%set(c, 'YTick', 0:log10(tmax), 'YTickLabel', sprintfc('10^{%d}',0:log10(tmax)) );
set(gca, 'FontSize', 24, 'XTick', 1:4, 'YTick', []);
xlabel('Number of interactions');
%title('Maximum equilibration time');
title(sprintf('Fraction with t_{eq} > %d', t_min));
%ylabel(c, 'Time steps', 'Interpreter', 'tex')
ylabel(c, 'Fraction of simulations', 'Interpreter', 'tex')

% Save figures
set(h, 'Units', 'Inches', 'Position', [0.1 0.1 12 8]);
qsave = 1;
%fname_str = sprintf('graph_plot_t_out_filtered_cond3_t_min_%d_jet', t_min);
fname_str = sprintf('graph_plot_t_out_filtered_cond3_t_min_%d_viridis', t_min);
save_figure(h, 12, 8, fullfile(save_folder, fname_str),'.pdf', qsave);
%% (3) Oscillation prevalence

% find x and y positions, order by value of output variable
x_graph = num_int;
y_graph = zeros(n_networks, 1);
output_var = frac_osc;
for ii=1:4
    idx = find(num_int==ii);
    output_var_temp = output_var(idx);
    [output_var_sorted, sort_idx] = sort(output_var_temp, 'ascend');
    y_graph(idx(sort_idx)) = linspace(-1, 1, numel(idx) );
end

% Plot oscillation prevalence
h3=figure;
hold on
G = graph(graph_edges);
% Layouts: circle, force(3), layered, subspace(3)
p = plot(G, 'Layout', 'force', 'MarkerSize', 1, 'LineWidth', 1.5, 'NodeCData',...
    frac_osc, 'NodeLabel', []);
% Plot marker edges for classes
scatter(x_graph(class{1}), y_graph(class{1}), 20^2, output_var(class{1}),...
    'o', 'filled');
scatter(x_graph(class{2}), y_graph(class{2}), 20^2, output_var(class{2}),...
    '^', 'filled');
scatter(x_graph(class{3}), y_graph(class{3}), 20^2, output_var(class{3}),...
    's', 'filled');

p.XData = x_graph;
p.YData = y_graph;
for ii=1:n_networks
   text(p.XData(ii)+0.05,p.YData(ii)-0.05, p.ZData(ii), num2str(ii), 'fontsize', 20, 'Color', 'k');
end
colormap(viridis);
%colormap(jet);

c=colorbar;
c.FontSize = 24;
caxis([0 1]);

set(gca, 'FontSize', 24, 'XTick', 1:4, 'YTick', [], 'ZTick', []);
xlabel('Number of interactions');
title('Oscillation prevalence')

% Save figure
qsave = 1;
fname_str = 'graph_plot_osc_prevalence_rearranged_marked_v4_viridis';
save_figure(h3, 12, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

%% (4) Final spatial order
% find x and y positions, order by value of output variable
x_graph = num_int;
y_graph = zeros(n_networks, 1);
output_var = frac_spatially_ordered;
for ii=1:4
    idx = find(num_int==ii);
    output_var_temp = output_var(idx);
    [output_var_sorted, sort_idx] = sort(output_var_temp, 'ascend');
    y_graph(idx(sort_idx)) = linspace(-1, 1, numel(idx) );
end

% Plot graph
h2 = figure;
hold on
G = graph(graph_edges);
% Layouts: circle, force(3), layered, subspace(3)
p = plot(G, 'Layout', 'force', 'MarkerSize', 1, 'LineWidth', 1.5, 'NodeCData',...
    zeros(n_networks, 1), 'NodeLabel', []);
% Plot marker edges for classes
scatter(x_graph(class{1}), y_graph(class{1}), 20^2, frac_spatially_ordered(class{1}),...
    'o', 'filled');
scatter(x_graph(class{2}), y_graph(class{2}), 20^2, frac_spatially_ordered(class{2}),...
    '^', 'filled');
scatter(x_graph(class{3}), y_graph(class{3}), 20^2, frac_spatially_ordered(class{3}),...
    's', 'filled');
p.XData = x_graph;
p.YData = y_graph;
for ii=1:n_networks
	text(p.XData(ii)+0.05,p.YData(ii)-0.05, p.ZData(ii), num2str(ii), 'fontsize', 20, 'Color', 'k');
end
%colormap(jet);
colormap(viridis);
c=colorbar;
c.FontSize = 24;

caxis([0 0.1]);
set( c, 'YTick', 0:0.05:0.1 );
set(gca, 'FontSize', 24, 'XTick', 1:4, 'YTick', []);
xlabel('Number of interactions');
title('Fraction spatially ordered')

qsave = 1;
fname_str = strrep(sprintf('graph_plot_frac_I_final_geq_%.1f_mean_rearranged_marked_v4_viridis', I_min), '.', 'p');
save_figure(h2, 12, 8, fullfile(save_folder, fname_str),'.pdf', qsave);
%% Correlation between final spatial order and oscillation prevalence
h = figure;
hold on
%
% Plot for different "classes" -> no clear trend
p = cell(3,1);
p{1} = scatter(frac_osc(class{1}), frac_spatially_ordered(class{1}), 100, 'filled', 'ko');
p{2} = scatter(frac_osc(class{2}), frac_spatially_ordered(class{2}), 100, 'filled', 'k^');
p{3} = scatter(frac_osc(class{3}), frac_spatially_ordered(class{3}), 100, 'filled', 'ks');
legend([p{1} p{2} p{3}], {'Complex', 'Oscillatory', 'Non-periodic'}, 'Location', 'nw');
%}
%{
% Plot against number of interactions -> no trend
scatter(frac_osc(x_graph==1), frac_spatially_ordered(x_graph==1), 100, 'filled', 'o');
scatter(frac_osc(x_graph==2), frac_spatially_ordered(x_graph==2), 100, 'filled', 'o');
scatter(frac_osc(x_graph==3), frac_spatially_ordered(x_graph==3), 100, 'filled', 'o');
scatter(frac_osc(x_graph==4), frac_spatially_ordered(x_graph==4), 100, 'filled', 'o');
%}
% Plot against number of repressive interactions 
%{
p = cell(5, 1);
cmap = get(gca, 'colororder'); %viridis(5);
for ii=0:4 %4
    % total # repressive interactions
    p{ii+1} = scatter(frac_osc(num_repressors==ii),...
        frac_spatially_ordered(num_repressors==ii), 100, cmap(ii+1,:), 'filled', 'o');
    % # repressive self-interactions
    %p{ii+1} = scatter(frac_osc(num_repr_self==ii),...
    %    frac_spatially_ordered(num_repr_self==ii), 100, cmap(ii+1,:), 'filled', 'o');
end
leg = legend([p{1} p{2} p{3}], sprintfc('%d', 0:2), 'Location', 'nw');
%title(leg,'# self-repressions')
title(leg,'# repressions')
leg.Title.Visible = 'on';
%legend([p{1} p{2} p{3} p{4} p{5}], sprintfc('%d', 0:4), 'Location', 'nw');
%}

% Plot averages over number of repressors
%{
frac_osc_by_numR_mean = zeros(5, 1);
frac_osc_by_numR_std = zeros(5, 1);
frac_spatially_ordered_by_numR_mean = zeros(5, 1);
frac_spatially_ordered_by_numR_std = zeros(5, 1);

network_counts = zeros(5,1);
cmap = get(gca, 'colororder'); %viridis(5);
cmap = cmap(1:5, :);
legend_str = cell(5, 1);
for ii=0:4
    p{ii+1} = scatter( mean(frac_osc(num_repressors==ii)), ...
        mean(frac_spatially_ordered(num_repressors==ii)),...
        500, cmap(ii+1, :), '^', 'filled');
    frac_osc_by_numR_mean(ii+1) = mean(frac_osc(num_repressors==ii));
    frac_osc_by_numR_std(ii+1) = std( frac_osc(num_repressors==ii), 1 );
    frac_spatially_ordered_by_numR_mean(ii+1) = mean(frac_spatially_ordered(num_repressors==ii));
    frac_spatially_ordered_by_numR_std(ii+1) = std( frac_spatially_ordered(num_repressors==ii), 1 );
    network_counts(ii+1) = sum(num_repressors==ii);
    legend_str{ii+1} = sprintf('n_{R}=%d (n=%d)', ii, network_counts(ii+1));
end
errorbar(frac_osc_by_numR_mean, frac_spatially_ordered_by_numR_mean, ...
    frac_spatially_ordered_by_numR_std, frac_spatially_ordered_by_numR_std, ... % y errorbars
    frac_osc_by_numR_std, frac_osc_by_numR_std); % x errorbars
%legend([p{1}, p{2}, p{3}, p{4}, p{5}], sprintfc('%d', 0:4), 'FontSize', 20, 'Location', 'nw');
legend([p{1}, p{2}, p{3}, p{4}, p{5}], legend_str, 'FontSize', 20, 'Location', 'nw');
%scatter(frac_osc_mean_by_numR, frac_spatially_ordered_mean_by_numR, 200, cmap, '^', 'filled');
%}

xlim([0 1]);
ylim([0 0.1]);
xlabel('Fraction periodic');
ylabel('Fraction spatially ordered');
set(gca, 'FontSize', 20, 'XTick', 0:0.2:1, 'YTick', 0:0.02:0.1);
%title(sprintf('$I>%.1f$', I_min ))

% label networks
%
for ii=1:n_networks
    % gate data
    cond = frac_osc(ii)>0.1 && frac_spatially_ordered(ii)>0.007;
    if cond || ismember(ii, class{1})
        text(frac_osc(ii)+0.01, frac_spatially_ordered(ii)-0.002, num2str(ii), 'fontsize', 16, 'Color', 'k');
    end
end
%}

qsave = 0;
save_folder2 = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2';
% networks, by class
fname_str = strrep(sprintf('frac_osc_vs_frac_ordered_by_network_class_I_min_%.1f', I_min), '.', 'p');
% networks, by # repressors
%fname_str = strrep(sprintf('frac_osc_vs_frac_ordered_by_network_num_repressors_I_min_%.1f', I_min), '.', 'p');
% averages, by # repressors
%fname_str = strrep(sprintf('frac_osc_vs_frac_ordered_by_num_repressors_mean_I_min_%.1f_errorbars', I_min), '.', 'p');
% averages, by # repressive self-interactions
%fname_str = strrep(sprintf('frac_osc_vs_frac_ordered_by_num_repressive_self_int_I_min_%.1f', I_min), '.', 'p');
save_figure(h, 10, 8, fullfile(save_folder2, fname_str),'.pdf', qsave);
%% Correlation part 2: testing hypotheses
% Is the positive correlation due to a direct relation between periodicity
% and spatial order? 

% Reload trajectories, get indices of spatially ordered and non-ordered
% trajectories
I_min = 0.3;
networks_sel = 1:44;

% output
count_periodic = zeros( numel(networks_sel), 1 ); % #periodic trajectories
count_non_periodic = zeros( numel(networks_sel), 1 ); % #non-periodic trajectories, not reaching tmax
num_ordered_p = zeros( numel(networks_sel), 1 ); % #ordered, periodic trajectories
num_ordered_non_p = zeros( numel(networks_sel), 1 ); % #ordered, non-periodic trajectories
% filter by period
% periods_chosen = [2, 3, 4, Inf, multiple of gz, other]
count_by_period = zeros( numel(networks_sel), 6 ); % #spatially ordered trajectories, by period
num_by_period = zeros( numel(networks_sel), 6 ); % total # number of trajectories, by period

for network=networks_sel
    % filter out trajectories that reached tmax
    %disp(sum(sum(t_out_all_net==tmax, 3), 2)); % # trajectories reached tmax per network
    idx_tmax_reached = t_out_all_net(network,:,:)==tmax;
    
    % check periodicity of ordered and disordered trajectories
    this_periods = period_all_net(network, :, :);
    periodic_idx = (this_periods<Inf);
    count_periodic(network) = sum(periodic_idx(:));
    count_non_periodic(network) = sum(sum(~periodic_idx & ~idx_tmax_reached));
    
    % Get final I of networks
    %network = 1;
    load_folder = 'H:\My Documents\Multicellular automaton\data\two_signals\batch_sim_all_topologies\analysis_initial_p_I';
    fname_str = sprintf('batch_sim_analyzed_data_batch2_p_I_final_network_%d_old%d',...
        network, network_all(network) );
    fname = fullfile(load_folder, fname_str);
    disp(fname);
    load(fname, 'I_final_network');
    
    % Separate into I_final for periodic and non-periodic trajectories
    I_final_network = reshape(I_final_network, size(I_final_network, 1)*size(I_final_network, 2), 2);
    %I_final_periodic = I_final_network(periodic_idx, :);
    %I_final_non_periodic = I_final_network(~periodic_idx & ~idx_tmax_reached, :);
    
    % Find fraction of spatially ordered among each of the two 
    num_ordered_p(network) = sum(...
        I_final_network(periodic_idx,1)>I_min &...
        I_final_network(periodic_idx,2)>I_min);
    num_ordered_non_p(network) = sum(...
        I_final_network(~periodic_idx & ~idx_tmax_reached, 1)>I_min &...
        I_final_network(~periodic_idx & ~idx_tmax_reached, 2)>I_min);
    
    % #spatially ordered trajectories per period  
    % clear periodic_idx   % clear variables for memory
    idx_I = (I_final_network(:,1) > I_min & I_final_network(:,2) > I_min & ~idx_tmax_reached(:));
    periods_ordered = this_periods(idx_I);
    count_by_period_temp = histcounts( categorical(periods_ordered), {'2', '3', '4', 'Inf'} );
    count_by_period_temp(5) = sum( mod(periods_ordered, gz)==0 );
    count_by_period_temp(6) = numel(periods_ordered) - sum(count_by_period_temp(1:5));
    count_by_period(network, :) = count_by_period_temp;
    
    % Total # number of trajectories, by period 
    num_by_period_temp = histcounts( categorical(this_periods), {'2', '3', '4', 'Inf'} );
    num_by_period_temp(5) = sum( mod(this_periods(:), gz)==0 );
    num_by_period_temp(6) = numel(this_periods) - sum(num_by_period_temp(1:5));
    num_by_period(network, :) = num_by_period_temp;    
    
    % other method (reverse): check periods of spatially ordered trajectories
    %{
    % get indices of ordered trajectories
    this_indices = (I_final_network(:,:,1)>I_min & I_final_network(:,:,2)>I_min);
    
    %this_periods = period_all_net(network, :, :);
    %temp = this_periods(this_indices);
    %temp2 = this_periods(~this_indices);
    %}
end
frac_ordered_p = num_ordered_p./count_periodic;
frac_ordered_non_p = num_ordered_non_p./(count_non_periodic);


%% Plot fraction spatially ordered, averaged over periodic/non-periodic only
% Ordered by network
h=figure;
hold on
barwidth = 1;
bar(1:n_networks, [frac_ordered_p frac_ordered_non_p], barwidth); %, 'bo-');
legend({'X = Periodic', 'X = Non-periodic'});
%bar(40:44, [frac_ordered_p(40:44) frac_ordered_non_p(40:44)], 1); %, 'bo-');
%barh(1:n_networks, ); %, 'ro-');
set(gca, 'XTick', 1:n_networks);
set(gca, 'FontSize', 20);
xlabel('Network');
%ylabel('Fraction spatially ordered');
ylabel('$P($ordered$|$X$)$', 'Interpreter', 'latex');

qsave = 0;
save_folder2 = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2';
fname_str = strrep(sprintf('frac_ordered_network_bar_I_min_%.1f', I_min), '.', 'p');
save_figure(h, 24, 6, fullfile(save_folder2, fname_str),'.pdf', qsave);

% Ordered by frequency
frac_ordered_p(isnan(frac_ordered_p)) = 0;
[~, sort_idx] = sort(frac_ordered_p, 'ascend');
h=figure;
bar( 1:n_networks, [frac_ordered_p(sort_idx) frac_ordered_non_p(sort_idx)]);
legend({'X = Periodic', 'X = Non-periodic'}, 'Location', 'nw');
set(gca, 'XTick', 1:n_networks, 'XTickLabel', sort_idx );
set(gca, 'FontSize', 20);
xlabel('Network');
%ylabel('Fraction spatially ordered');
ylabel('$P($ordered$|$X$)$', 'Interpreter', 'latex');

qsave = 0;
save_folder2 = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2';
fname_str = strrep(sprintf('frac_ordered_network_bar_I_min_%.1f_v2_sorted', I_min), '.', 'p');
save_figure(h, 24, 6, fullfile(save_folder2, fname_str),'.pdf', qsave);

%% Plot mean|Complex and mean|Oscillatory per class 
%{
h = figure;
hold on
bar([1 2], [mean(frac_ordered_p(class{1})) mean(frac_ordered_p(class{2}))]);
errorbar( 1, mean(frac_ordered_p(class{1})), std(frac_ordered_p(class{1})) );
errorbar( 2, mean(frac_ordered_p(class{2})), std(frac_ordered_p(class{2})) );
set(gca, 'XTick', [1 2], 'XTickLabels', {'Complex', 'Oscillatory'});
% -> Not statistically significant
%}
%% Plot fraction spatially ordered, averaged over all trajectories
% Sort by frequency
frac1 = num_ordered_p./(count_periodic+count_non_periodic);
frac2 = num_ordered_non_p./(count_periodic+count_non_periodic);
frac1(isnan(frac1)) = 0;
[~, sort_idx] = sort(frac1+frac2, 'ascend');

h = figure;
barwidth = 1;
bar(1:n_networks, [frac1(sort_idx) frac2(sort_idx)], barwidth, 'stacked'); %, 'bo-');
legend({'Periodic', 'Non-periodic'}, 'Location', 'nw');
%bar(40:44, [frac_ordered_p(40:44) frac_ordered_non_p(40:44)], 1); %, 'bo-');
%barh(1:n_networks, ); %, 'ro-');
set(gca, 'XTick', 1:n_networks,  'XTickLabel', sort_idx );
set(gca, 'FontSize', 20);
xlabel('Network');
ylabel('Fraction spatially ordered');

qsave = 0;
save_folder2 = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2';
fname_str = strrep(sprintf('frac_ordered_network_bar_I_min_%.1f_v3', I_min), '.', 'p');
save_figure(h, 24, 6, fullfile(save_folder2, fname_str),'.pdf', qsave);

%% Analyze spatially ordered periodic trajectories by period

% Which periods most likely lead to spatial order?
frac_by_period = count_by_period./num_by_period; % fraction of spatially ordered trajectories per period
idx_nan = isnan(frac_by_period);
frac_by_period(idx_nan) = 0;
%
% sort data
sort_idx_final = 1:n_networks;
for ii=[4 1 2 3 6 5] %sort data one by one
    [~, sort_idx] = sort( frac_by_period(:,ii), 'ascend' );
    frac_by_period = frac_by_period(sort_idx, :);
    sort_idx_final = sort_idx_final(sort_idx);
end
%}
frac_by_period(idx_nan(sort_idx_final,:)) = NaN; % set the non-existing fractions to NaN again
h=figure;
hold on
plot(frac_by_period, 'o', 'LineWidth', 2, 'MarkerSize', 10);
%scatter(repmat((1:n_networks)', 1, 6), frac_by_period, 'filled');
legend({'2', '3', '4', 'Inf', 'mult. gz', 'other'}, 'Location', 'nw');
xlabel('Network');
%ylabel('Fraction spatially ordered');
ylabel('$P($ordered$|$period$)$', 'Interpreter', 'latex');
set(gca, 'XTick', 1:n_networks, 'XTickLabel', sort_idx_final );
set(gca, 'YTick', 0:0.2:1);
set(gca, 'FontSize', 20);
%set(h, 'Units', 'Inches', 'Position', [1 1 24 6]);
box on

qsave = 0;
save_folder2 = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2';
fname_str = strrep(sprintf('frac_ordered_by_period_network_I_min_%.1f_plot_v2', I_min), '.', 'p');
save_figure(h, 24, 6, fullfile(save_folder2, fname_str),'.pdf', qsave);

% How do the spatially ordered trajectories break down by period?

%% Select complex class networks only, examine effects
frac_period_above_4 = sum(count_by_period(class{1}, 5:6), 2)./sum(num_by_period(class{1}, 5:6), 2);
frac_period_mult_gz = count_by_period(class{1}, 5)./num_by_period(class{1}, 5);

h=figure;
hold on
bar(1:numel(class{1}), [frac_period_mult_gz frac_period_above_4]);
P_max = max(frac_spatially_ordered); % maximumm fraction of spatially ordered trajectories, over all networks and trajectories
plot([0 numel(class{1})+1], [P_max P_max], 'k--');
legend({'Mult. gz', 'Period>4'}, 'Location', 'se');
xlabel('Network');
%ylabel('Fraction spatially ordered');
ylabel('$P($ordered$|$period$)$', 'Interpreter', 'latex');
set(gca, 'XTick', 1:numel(class{1}), 'XTickLabel', class{1});
set(gca, 'FontSize', 20);
%set(h, 'Units', 'Inches', 'Position', [1 1 24 6]);

for ii=1:numel(class{1})
    % Label number of period mult. gz trajectories
    text( ii-0.3, frac_period_mult_gz(ii) + 0.02,...
        sprintf('%d', count_by_period(class{1}(ii), 5)), 'fontsize', 16);
    % Label number of period mult. gz trajectories
    text( ii, frac_period_above_4(ii) + 0.02,...
        sprintf('%d', sum(count_by_period(class{1}(ii), 5:6), 2) ) , 'fontsize', 16);
end
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);

qsave = 0;
save_folder2 = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2';
fname_str = strrep(sprintf('frac_ordered_by_period_complex_network_I_min_%.1f_bar_v2_with_counts', I_min), '.', 'p');
save_figure(h, 10, 8, fullfile(save_folder2, fname_str),'.pdf', qsave);


%% Plot periods as heat map (not as graph)
uniq_periods = unique(period_all_net);

% counts = zeros( n_networks, numel(uniq_periods) );
edges = [2:4 5:5:100 Inf];
nbins = numel(edges)-1;
counts = zeros( n_networks, numel(edges)-1 );
for ii=1:n_networks
    %C = categorical(period_all_net(ii,:), uniq_periods);
    %counts(ii,:) = histcounts(C);
    
    counts(ii,:) = histcounts(period_all_net(ii,:), edges);
end

h=figure;
imagesc([log10(counts(class{1}, 1:end));...
    log10(counts(class{2}, 1:end));...
    log10(counts(class{3}, 1:end))])

% xlabels
xlabels = cell(nbins, 1);
xlabels(1:3) = sprintfc('%d', edges(1:3));
for ii=4:nbins
    xlabels{ii} = sprintf('%d-%d', edges(ii), edges(ii+1));
end
set(gca, 'XTick', 1:nbins, 'XTickLabel', xlabels);

% ylabels
set(gca, 'YTick', 1:44, 'YTickLabel', [class{1} class{2} class{3}]);

% Save figure
qsave = 0;
fname_str = 'periods_density_by_network_imagesc';
save_figure(h, 15, 10, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Histograms by class
uniq_periods = unique(period_all_net);

counts = zeros( 3, numel(uniq_periods) );
edges = [2:4 5:5:100 Inf];
nbins = numel(edges)-1;
counts = zeros( 3, numel(edges)-1 );
for ii=1:3
    %C = categorical(period_all_net(class{ii},:), uniq_periods);
    %counts(ii,:) = histcounts(C);
        
    counts(ii,:) = histcounts(period_all_net(ii,:), edges);
end

h=figure;
hold on
bar3( log10(counts) );
%bar3(1*ones(nbins ,1), log10( counts(1,:) ))
%bar3(2*ones(nbins ,1), log10( counts(2,:) ))
%bar3(3*ones(nbins ,1), log10( counts(3,:) ))

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