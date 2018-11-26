%% Obtain statistics on all possible topologies by simulation, plot results as graphs
clear variables
close all
clc
set(0,'defaulttextinterpreter','latex');
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
class{1} = [15	16	19	20	32	33	34	36	43]; % complex dynamics
class{2} = [2	5	8	10	11	12	13	14	18	21	22	23	25	27, ...
                29	30	31	35	37	38	39	40	41	42	44]; % simple oscillations
class{3} = [1	3	4	6	7	9	17	24	26	28]; % no oscillations

% Networks by number of repressive interactions
num_repressors = zeros(n_networks, 1);
for ii=1:n_networks
    num_repressors(ii) = sum(M_int_all_reduced{ii}(:)==-1);
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
%% (1) Plot max(t_out)
t_out_max = max(t_out_all_net(:,:), [], 2);
t_out_max_log10 = log10(t_out_max);

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
   text(p.XData(ii)+0.06,p.YData(ii)-0.05, p.ZData(ii), num2str(ii), 'fontsize', 20, 'Color', 'k');
end
%colormap(viridis);
colormap(jet);
c=colorbar;
c.FontSize = 24;
caxis([0 log10(tmax)]);
set(c, 'YTick', 0:log10(tmax), 'YTickLabel', sprintfc('10^{%d}',0:log10(tmax)) );
set(gca, 'FontSize', 24, 'XTick', 1:4, 'YTick', []);
xlabel('Number of interactions');
ylabel(c, 'max($$t_{eq}$$)', 'Interpreter', 'latex')

% Save figures
set(h1, 'Units', 'Inches', 'Position', [0.1 0.1 12 8]);
qsave = 0;
fname_str = 'graph_plot_t_out_max_log_rearranged_marked_v3_jet_size_12_8';
save_figure(h1, 12, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

%% (2) Plot t_out mean
t_out_mean = mean(t_out_all_net(:,:), 2);

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
	text(p.XData(ii)+0.05,p.YData(ii)-0.05, p.ZData(ii), num2str(ii), 'fontsize', 16, 'Color', 'k');
end
%colormap(viridis);
colormap(jet);
c=colorbar;
c.FontSize = 20;
%caxis([0 250]);
caxis([0 3]);
set(c, 'YTick', 0:3, 'YTickLabel', sprintfc('10^{%d}',0:3) );
set(gca, 'FontSize', 20, 'XTick', 1:4, 'YTick', []);
xlabel('Number of interactions');
title('$$\langle t_{out} \rangle$$')

qsave = 1;
fname_str = 'graph_plot_t_out_mean_rearranged_marked_v2_jet';
save_figure(h2, 15, 10, fullfile(save_folder, fname_str),'.pdf', qsave);

%% (3) Oscillation prevalence
frac_osc = zeros(n_networks, 1); % Fraction of parameter sets which are capable of generating oscillations
frac_osc_2 = zeros(n_networks, 1); % Net fraction of simulations which generate oscillations
for ii=1:n_networks
    period_temp = squeeze(period_all_net(ii,:,:));
    [idx1, idx2] = find(period_temp~=Inf);
    frac_osc(ii) = numel(unique(idx1))/n_pset;
    frac_osc_2(ii) = numel(idx1)/(n_pset*nsim);
end

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
   text(p.XData(ii)+0.05,p.YData(ii)-0.05, p.ZData(ii), num2str(ii), 'fontsize', 16, 'Color', 'k');
end
%colormap(viridis);
colormap(jet);

c=colorbar;
c.FontSize = 20;
caxis([0 1]);

set(gca, 'FontSize', 20, 'XTick', 1:4, 'YTick', [], 'ZTick', []);
xlabel('Number of interactions');
title('Oscillation prevalence')

% Save figure
qsave = 0;
fname_str = 'graph_plot_osc_prevalence_rearranged_marked_v2_jet';
save_figure(h3, 15, 10, fullfile(save_folder, fname_str),'.pdf', qsave);

%% (4) Final spatial order
I_min = 0.3;

% Load files, find fraction of networks giving spatial order
networks_sel = 1:44;
num_ordered = zeros( numel(networks_sel), 1 );
for network=networks_sel
    %network = 1;
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
	text(p.XData(ii)+0.05,p.YData(ii)-0.05, p.ZData(ii), num2str(ii), 'fontsize', 16, 'Color', 'k');
end
colormap(jet);
c=colorbar;
c.FontSize = 20;

caxis([0 0.1]);
set(c, 'YTick', 0:0.05:0.1 );
set(gca, 'FontSize', 20, 'XTick', 1:4, 'YTick', []);
xlabel('Number of interactions');
title('Fraction spatially ordered')

qsave = 0;
fname_str = strrep(sprintf('graph_plot_frac_I_final_geq_%.1f_mean_rearranged_marked', I_min), '.', 'p');
save_figure(h2, 15, 10, fullfile(save_folder, fname_str),'.pdf', qsave);
%% Correlation between final spatial order and oscillation prevalence
h = figure;
hold on
%{
% Plot for different "classes" -> no clear trend
scatter(frac_osc(class{1}), frac_spatially_ordered(class{1}), 100, 'filled', 'ko');
scatter(frac_osc(class{2}), frac_spatially_ordered(class{2}), 100, 'filled', 'k^');
scatter(frac_osc(class{3}), frac_spatially_ordered(class{3}), 100, 'filled', 'ks');
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
for ii=0:4
    p{ii+1} = scatter(frac_osc(num_repressors==ii),...
        frac_spatially_ordered(num_repressors==ii), 100, cmap(ii+1,:), 'filled', 'o');
end
legend([p{1} p{2} p{3} p{4} p{5}], sprintfc('%d', 0:4), 'Location', 'nw');
%}

% Plot averages over number of repressors
%
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
%{
for ii=1:n_networks
    % gate data
    if frac_osc(ii)>0.1 && frac_spatially_ordered(ii)>0.007
        text(frac_osc(ii)+0.01, frac_spatially_ordered(ii)-0.002, num2str(ii), 'fontsize', 16, 'Color', 'k');
    end
end
%}

qsave = 1;
save_folder2 = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2';
% networks, by class
%fname_str = strrep(sprintf('frac_osc_vs_frac_ordered_by_network_class_I_min_%.1f', I_min), '.', 'p');
% networks, by # repressors
%fname_str = strrep(sprintf('frac_osc_vs_frac_ordered_by_network_num_repressors_I_min_%.1f', I_min), '.', 'p');
% averages, by # repressors
fname_str = strrep(sprintf('frac_osc_vs_frac_ordered_by_num_repressors_mean_I_min_%.1f_errorbars', I_min), '.', 'p');
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
qsave = 1;
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