%% Analyze data from batch simulations with OR logic
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'tex');

%% Load analyzed data
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\batch_sim_all_topologies_OR_logic_run1\analyzed_data';
%save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\batch_sim_all_topologies_OR_logic_run1\analyzed_data_1sim';
fname_str = 'batch_sim_OR_logic_all_topologies_analyzed_all';
%fname_str = 'batch_sim_OR_logic_all_topologies_analyzed_network_3';
fname_out = fullfile(save_folder, fname_str);
load(fname_out);

%% Basic parameters
n_pset = 10^4;
nsim = 10;
tmax = 10^4; 

save_folder_fig = fullfile('H:\My Documents\Multicellular automaton\figures', ...
    'batch_sim_all_topologies_OR_logic_run1');

%% Calculate for network graph representation 
% Calculate network distances (Hamming distance)
n_networks = numel(network_all);
graph_dist = zeros(n_networks);
for i1=1:n_networks
    for i2=1:n_networks
        %graph_dist(i,j) = 1;
        n_perm = 2;
        %M_sim_v1 = zeros(2); % measures how similar the two networks are
        M_dist = zeros(2); % measures how similar the two networks are
        for idx1=1:n_perm
            for idx2=1:n_perm
                M_int_1 = M_perm(M_int_all_reduced{i1}, idx1);
                M_int_2 = M_perm(M_int_all_reduced{i2}, idx2);
                %M_sim_v1(idx1, idx2) = sum(M_int_1(:)==M_int_2(:));
                M_dist(idx1, idx2) = sum(abs(M_int_1(:) - M_int_2(:)));
                %disp(M_int_1);
                %disp(M_int_2);
                %disp(M_sim(idx1, idx2));
            end
        end
        %graph_sim(i,j) = max(M_sim_v1(:));
        graph_dist(i1,i2) = min(M_dist(:));
    end
end
%graph_edges = (graph_sim==3); % If three links are the same, the networks are neighbours
graph_edges = (graph_dist==1); % If three links are the same, the networks are neighbours

% Calculate x, y positions of networks
% get the number of interactions of each network
num_int = zeros(n_networks, 1); 
for i1=1:n_networks
    num_int(i1) = sum(sum( abs(M_int_all_reduced{i1}) ));
end

% number of networks with a given # of interactions
num_net_by_int = zeros(4, 1);
for i1=1:4
    num_net_by_int(i1) = sum(num_int == i1);
end

%% Calculate data for later usage
% tmax
t_out_max = max(t_out_all_network(:,:), [], 2);
t_out_max_log10 = log10(t_out_max);

% t_mean
t_out_mean = mean(t_out_all_network(:,:), 2);

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
t_out_filtered = mean(t_out_all_network(:,:)>t_min, 2);

% Fraction oscillators
frac_osc = zeros(n_networks, 1); % Fraction of parameter sets which are capable of generating oscillations
frac_osc_2 = zeros(n_networks, 1); % Net fraction of simulations which generate oscillations
%ii = 15;
for ii=1:n_networks
    period_temp = squeeze(period_all_network(ii,:,:));
    [idx1, ~] = find(period_temp~=Inf);
    frac_osc(ii) = numel(unique(idx1))/n_pset; % <- check whether this is correct!
    frac_osc_2(ii) = numel(idx1)/(n_pset*nsim);
end

% Fraction spatially ordered
I_min = 0.3;
networks_sel = 1:44;
num_ordered = zeros( numel(networks_sel), 1 );
for network=networks_sel
    I_final_network = squeeze(I_final_all_network(network, :, :, :));
    %num_ordered(network) = sum(sum(I_final_network(:,:,1)>I_min & I_final_network(:,:,2)>I_min));
    %--------------TEMP correction for n_sim=1 (revert back later)---------
    num_ordered(network) = sum(sum(I_final_network(:,1)>I_min & I_final_network(:,2)>I_min));
    %----------------------------------------------------------------------
end
frac_spatially_ordered = zeros(n_networks, 1);
frac_spatially_ordered(networks_sel) = num_ordered(networks_sel)/(n_pset*nsim);

% Network classification
% Based on following results, propose classification
class = cell(3, 1);
class{1} = [15	16	19	20	33	34	36	43]; % complex dynamics
class{2} = [2	5	8	10	11	12	13	14	18	21	22	23	25	27, ...
                29	30	31	32  35	37	38	39	40	41	42	44]; % simple oscillations
class{3} = [1	3	4	6	7	9	17	24	26	28]; % no oscillations
%% -----------------Topologies capable of producing oscillations-----------
[idx1, idx2] = find(period_all_network(:,:)~=Inf);
n_osc = numel(unique(idx1)); % number of topologies capable of producing oscillations

h = figure;
plotHandle = pie([n_osc n_networks-n_osc]);
hText = findobj(plotHandle, 'Type', 'text'); % text object handles
percentValues = get(hText, 'String'); % percent values
countValues = string([n_osc; n_networks-n_osc]);
txt = {'Oscillatory: '; 'Non-oscillatory: '};
combinedtxt = strcat(txt, percentValues);
combinedtxt2 = strcat( '(', countValues, ' networks)');
oldExtents_cell = get(hText, 'Extent'); 
oldExtents = cell2mat(oldExtents_cell); % numeric array
hText(1).String = {combinedtxt{1}, combinedtxt2{1}};
hText(2).String = {combinedtxt{2}, combinedtxt2{2}};
set(plotHandle(2:2:end), 'FontSize', 20);

newExtents_cell = get(hText, 'Extent'); % cell array
newExtents = cell2mat(newExtents_cell); % numeric array 
width_change = (newExtents(:,3)-oldExtents(:,3));
signValues = sign(oldExtents(:,1));
offset = 0; %signValues.*(width_change/2);
textPositions_cell = get(hText,{'Position'}); % cell array
textPositions = cell2mat(textPositions_cell); % numeric array
textPositions(:,1) = textPositions(:,1) + offset; % add offset 

hText(1).Position = textPositions(1,:);
hText(2).Position = textPositions(2,:);

qsave = 0;
fname_str = 'oscillation_fraction_pie';
save_figure(h, 10, 8, fullfile(save_folder_fig, fname_str),'.pdf', qsave);

% Get list of topologies producing oscillations
topologies_osc = unique(idx1);
topologies_non_osc = setdiff(1:n_networks, unique(idx1));

%% Topologies capable of producing complex periods
% selected networks to include
% default
%network_sel = 1:numel(network_all);
% subset
% network_sel = class{1}; %reduced index
% custom
%
network_sel = [15 19 36 33 34 16 20 43]; % [15 16 19 20 33 34 36 43]; %reduced index
%period_all_net = period_all(network_sel_orig, :, :);
%}

h=figure;
hold on
box on
% sort networks by number of different periods (works only if all networks
% included)
%{
period_count = zeros(numel(network_sel), 1);
for i=1:numel(network_sel)
    period_count(i) = numel(unique(period_all_network(i,:)));
end
[~, idx_sort] = sort(period_count, 'descend');
for i=1:numel(network_sel)
    scatter( unique(period_all_network(idx_sort(i),:)), ...
        repmat(i, period_count(idx_sort(i)), 1), 35, 'kd', 'filled' );
end
set(gca, 'YTick', 1:numel(network_sel), 'YTickLabels', network_sel(idx_sort));
%}
% don't sort networks
%
for i=1:numel(network_sel)
    network_idx = network_sel(i);
    scatter( unique(period_all_network(network_idx,:)), ...
        repmat(i, numel(unique(period_all_network(network_idx,:))), 1), 35, 'kd', 'filled' );
end
set(gca, 'YTick', 1:numel(network_sel), 'YTickLabels', network_sel);
%}

% set other properties
set(gca, 'XScale', 'log');
xlim([1 10^4]);
ylim([0 numel(network_sel)+1]);
colors = parula(4);
plot([2 2], [0 numel(network_sel)+1], '--', 'Color',  colors(1,:) );
plot([3 3], [0 numel(network_sel)+1], '--', 'Color', colors(2,:));
plot([4 4], [0 numel(network_sel)+1], '--', 'Color', colors(3,:));
set(gca, 'FontSize', 20);
ylabel('Network', 'FontSize', 24);
xlabel('Period', 'FontSize', 24);

% Same topologies as those with high tmax!
qsave = 0;
fname_str = 'OR_logic_oscillation_periods_all_scatter_v2_subset';
save_figure(h, 10, 5, fullfile(save_folder_fig, fname_str),'.pdf', qsave);
%fname_str = 'OR_logic_oscillation_periods_all_scatter_v1_all_sorted';
%save_figure(h, 10, 8, fullfile(save_folder_fig, fname_str),'.pdf', qsave);

%% Distribution of periods among classes (Fig. S1C)
% Subdivide periodic trajectories into classes 
counts = cell(3,1);
cats = cell(3,1);
for i=1:3
    period_temp = squeeze(period_all_network(class{i}, :, :));
    n_network_temp = numel(class{i});
    n_osc_psets = 0;
    for j=1:n_network_temp
        [idx1, idx2] = find(squeeze(period_temp(j,:,:)~=Inf));
        n_osc_psets = n_osc_psets + numel(unique(idx1));
    end
    
    C = categorical(period_temp);
    [counts{i}, cats{i}] = histcounts(C);
    
    frac_osc(i) = n_osc_psets/(n_pset*n_network_temp); % count parameter sets
    %frac_osc_2(i) = numel(idx1)/(n_network_temp*n_pset*nsim); % count simulations
end

cats_new = cell(2,1);
cats_new{1} = {'2', '3', '4', '>4', 'Inf'};
counts_new = cell(2,1);
counts_new{1} = zeros(5, 1);
counts_new{1}(1) = counts{1}(1);
counts_new{1}(2) = counts{1}(2);
counts_new{1}(3) = counts{1}(3);
counts_new{1}(4) = sum(counts{1}(4:end-1));
counts_new{1}(5) = counts{1}(end);

cats_new{2} = cats{2}([1 2 3]); %{'2', '3', '4'};
counts_new{2} = zeros(5, 1);
counts_new{2}(1) = counts{2}(1);
counts_new{2}(3) = counts{2}(2);
counts_new{2}(5) = counts{2}(end);

% normalize by number of simulations
counts_normalized = cell(2,1);
counts_normalized{1} = counts_new{1}(1:end-1)/sum(counts_new{1}(1:end-1));
counts_normalized{2} = counts_new{2}(1:end-1)/sum(counts_new{2}(1:end-1));

%% Plot counts of new categories
% Dynamical spatial pattern producing networks
h1=figure;

plotHandle = pie(counts_normalized{1}, {'', '', '', ''}); %, cats_new(1:end-1));
set(plotHandle(2:2:end), 'FontSize', 20);
% set pie colors
rgbmatrix = parula(4);
for k = 1:numel(counts_new{1})-1
  % Create a color for this sector of the pie
  pieColorMap = rgbmatrix(k,:); 
  % Apply the colors we just generated to the pie chart.
  set(plotHandle(k*2-1), 'FaceColor', pieColorMap);
end
legend(cats_new{1}(1:4), 'location', 'eo', 'FontSize', 20);

qsave = 0;
fname_str = 'OR_logic_Oscillation_breakdown_complex_networks_pie';
save_figure(h1, 10, 8, fullfile(save_folder_fig, fname_str),'.pdf', qsave);

%%
% Dynamical temporal pattern producing networks
% --> only works if run immediately after previous section
h2=figure;
plotHandle = pie(counts_normalized{2}, {'', '', '', ''}); %, cats_new(1:end-1));
set(plotHandle(2:2:end), 'FontSize', 20);
% set pie colors
rgbmatrix = rgbmatrix([1 3], :);
for k = 1:2
  % Create a color for this sector of the pie
  pieColorMap = rgbmatrix(k,:); 
  % Apply the colors we just generated to the pie chart.
  set(plotHandle(k*2-1), 'FaceColor', pieColorMap);
end
%legend(cats_new{2}([1 3]), 'location', 'eo', 'FontSize', 20);
legend(cats_new{2}, 'location', 'eo', 'FontSize', 20);

qsave = 0;
fname_str = 'OR_logic_Oscillation_breakdown_osc_networks_pie';
save_figure(h2, 10, 8, fullfile(save_folder_fig, fname_str),'.pdf', qsave);
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
box on
hold on
G = graph(graph_edges);
% Layouts: circle, force(3), layered, subspace(3)
p = plot(G, 'Layout', 'force', 'MarkerSize', 0.1, 'LineWidth', 1.5, 'NodeCData',...
    log10(t_out_max), 'NodeLabel', []);

% Plot marker edges for classes
%scatter(x_graph, y_graph, 20^2, t_out_max_log10, 'filled');
%
scatter(x_graph(class{1}), y_graph(class{1}), 20^2, t_out_max_log10(class{1}),...
    'o', 'filled');
scatter(x_graph(class{2}), y_graph(class{2}), 20^2, t_out_max_log10(class{2}),...
    '^', 'filled');
scatter(x_graph(class{3}), y_graph(class{3}), 20^2, t_out_max_log10(class{3}),...
    's', 'filled');
%}
p.XData = x_graph;
p.YData = y_graph;
for ii=1:n_networks
   text(p.XData(ii)+0.06,p.YData(ii)-0.05, p.ZData(ii), num2str(ii),...
       'fontsize', 20, 'Color', 'k', 'interpreter', 'tex');
end
cmap = 'jet';
colormap(cmap);
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
fname_str = sprintf('OR_logic_graph_plot_t_out_max_log_rearranged_marked_%s', cmap);
save_figure(h1, 12, 8, fullfile(save_folder_fig, fname_str),'.pdf', qsave);

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
box on
hold on
G = graph(graph_edges);
% Layouts: circle, force(3), layered, subspace(3)
p = plot(G, 'Layout', 'force', 'MarkerSize', 1, 'LineWidth', 1.5, 'NodeCData',...
    t_out_mean, 'NodeLabel', []);
% Plot marker edges for classes
%scatter(x_graph, y_graph, 20^2, log10(output_var), 'filled');
%
scatter(x_graph(class{1}), y_graph(class{1}), 20^2, log10(output_var(class{1})),...
    'o', 'filled');
scatter(x_graph(class{2}), y_graph(class{2}), 20^2, log10(output_var(class{2})),...
    '^', 'filled');
scatter(x_graph(class{3}), y_graph(class{3}), 20^2, log10(output_var(class{3})),...
    's', 'filled');
%}
p.XData = x_graph;
p.YData = y_graph;
for ii=1:n_networks
	text(p.XData(ii)+0.05,p.YData(ii)-0.05, p.ZData(ii), num2str(ii),...
        'fontsize', 16, 'Color', 'k', 'Interpreter', 'tex');
end
cmap = 'jet';
colormap(cmap);
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
fname_str = strcat('OR_logic_graph_plot_t_out_mean_rearranged_marked_', cmap);
save_figure(h2, 12, 8, fullfile(save_folder_fig, fname_str),'.pdf', qsave);

%% (3) Oscillation prevalence
% find x and y positions, order by value of output variable
x_graph = num_int;
y_graph = zeros(n_networks, 1);
%output_var = frac_osc;
output_var = frac_osc_2;
for ii=1:4
    idx = find(num_int==ii);
    output_var_temp = output_var(idx);
    [output_var_sorted, sort_idx] = sort(output_var_temp, 'ascend');
    y_graph(idx(sort_idx)) = linspace(-1, 1, numel(idx) );
end

% Plot oscillation prevalence
h3 = figure;
box on
hold on
G = graph(graph_edges);
% Layouts: circle, force(3), layered, subspace(3)
p = plot(G, 'Layout', 'force', 'MarkerSize', 1, 'LineWidth', 1.5, 'NodeCData',...
    frac_osc, 'NodeLabel', []);
% Plot marker edges for classes
%scatter(x_graph, y_graph, 20^2, output_var, 'filled');
%
scatter(x_graph(class{1}), y_graph(class{1}), 20^2, output_var(class{1}),...
    'o', 'filled');
scatter(x_graph(class{2}), y_graph(class{2}), 20^2, output_var(class{2}),...
    '^', 'filled');
scatter(x_graph(class{3}), y_graph(class{3}), 20^2, output_var(class{3}),...
    's', 'filled');
%}

p.XData = x_graph;
p.YData = y_graph;
for ii=1:n_networks
   text(p.XData(ii)+0.05,p.YData(ii)-0.05, p.ZData(ii), num2str(ii), 'fontsize', 20, 'Color', 'k');
end
cmap = 'jet';
%colormap(viridis);
colormap(cmap);
c=colorbar;
c.FontSize = 24;
caxis([0 1]);
ylabel(c, 'Fraction');

set(gca, 'FontSize', 24, 'XTick', 1:4, 'YTick', [], 'ZTick', []);
xlabel('Number of interactions');
title('Oscillation prevalence')

% Save figure
qsave = 0;
fname_str = sprintf('OR_logic_graph_plot_osc_prevalence_rearranged_marked_%s', cmap);
save_figure(h3, 12, 8, fullfile(save_folder_fig, fname_str),'.pdf', qsave);

%% Plot average final I1, I2

% Pre-processing
gz = 15;
networks_sel = 1:44;
num_classes = 4; % number of different classes of trajectories to consider
num_ordered = zeros( numel(networks_sel), 1 );
I_mean_all = zeros( numel(networks_sel), num_classes, 2 );
I_std_all = zeros( numel(networks_sel), num_classes, 2 );
for network=networks_sel
    I_final_network = squeeze(I_final_all_network(network, :, :, :));
    
    % number of ordered networks
    num_ordered(network) = sum(sum(I_final_network(:,:,1)>I_min & I_final_network(:,:,2)>I_min));
    
    % final spatial order
    idx_class = {};
    
    % ordering 3
    idx_class{1} = squeeze(period_all_network(network, :, :) == Inf);
    idx_class{2} = squeeze(period_all_network(network, :, :) < 5);
    idx_class{3} = squeeze(period_all_network(network, :, :) < Inf & ...
        period_all_network(network, :, :) >= 5 & mod(period_all_network(network, :, :),gz)~=0);
    idx_class{4} = squeeze(period_all_network(network, :, :) < Inf & ...
        mod(period_all_network(network, :, :),gz)==0);
    
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

%% Plot avg. I for different classes, per network
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
legend({'T=\infty', 'T<5', 'T>=5, mod(T,gz)~=0', 'mod(T,gz)=0'}, 'Location', 'nw');
%legend({'T=\infty', 'T<5', '5<=T<gz', 'T>=gz'}, 'Location', 'nw');

xlabel('Network');
ylabel(sprintf('Final I^{(%d)}', mol), 'Interpreter', 'tex');
set(gca, 'XTick', 1:n_networks, 'XTickLabel', sort_idx_final );
%set(gca, 'YTick', 0:0.2:1);
set(gca, 'FontSize', 20);
set(h, 'Units', 'Inches', 'Position', [1 1 24 6]);
box on
ylim([-0.2 0.8]);

% Save plot
qsave = 0;
fname_str = sprintf('final_I%d_mean_by_network_errorbar_classification3', mol);
fname = fullfile(save_folder_fig, fname_str);
save_figure(h, 24, 6, fname, '.pdf', qsave);

%% 


%% Functions
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