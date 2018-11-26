%% Obtain statistics on all possible topologies by simulation
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
%parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies';
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2';

%% Load saved data
%
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\batch_sim_all_topologies';
fname_str = sprintf('batch_sim_analyzed_data_batch2.mat');
load(fullfile(folder, fname_str));
n_networks = numel(M_int_all_reduced);
%}
%% ------------ Analyze loaded data---------------------------
%
set(0, 'defaulttextinterpreter', 'latex');
save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\batch_sim_all_topologies_run2';

% divide into subsets
network_set1 = [15	16	19	20	32	33	34	36	43]; % complex dynamics
network_set2 = [2	5	8	10	11	12	13	14	18	21	22	23	25	27, ...
                29	30	31	35	37	38	39	40	41	42	44]; % simple oscillations
network_set3 = [1	3	4	6	7	9	17	24	26	28]; % no oscillations

% transform back to original indices
network_set1_orig = network_all(network_set1); % complex dynamics
network_set2_orig = network_all(network_set2);
network_set3_orig = network_all(network_set3);
network_sets = {network_set1_orig, network_set2_orig, network_set3_orig};

%{
network_conversion = zeros(numel(M_int_all), 1);
for i=1:numel(M_int_all)
    M_int = M_int_all{i};
    for j=1:numel(M_int_all_reduced)
        M_int_reduced = M_int_all_reduced{j};
        if all(M_int(:) == M_int_reduced(:))
            network_conversion(i) = j;
            break
        end
    end
end
%}
%% Maximum t_out per topology

% selected networks to include
% default
network_sel_orig = network_all;
network_sel = 1:numel(network_all);
% custom
%network_sel_orig = network_set3_orig; %original index
%network_sel = network_set3; %reduced index

t_out_all_net = t_out_all(network_sel_orig, :, :);

% Plot histogram of distribution of max(t_out)
%{
figure;
histogram( max(t_out_all_net(:,:), [], 2), 0:tmax/20:tmax)
%plot(1:n_networks, max(t_out_all_net(:,:), [], 2), 'bo');
xlabel('Maximum t_{out}');
ylabel('Number of network topologies');
%}

% Plot bar graph of max(t_out) for each topology
[t_out_max_sorted, sort_idx] = sort( max(t_out_all_net(:,:), [], 2), 'descend');
h=figure;
hold on
barh(t_out_max_sorted);
plot( [tmax tmax], [0 numel(network_sel)+1], 'r-');
plot( [1 1],[0 numel(network_sel)+1], 'r-');
set(gca, 'YTick', 1:numel(network_sel), 'YTickLabels', string(network_sel(sort_idx)));
set(gca, 'XScale', 'log');
xlim([10^(-0.1) 10^5]);
set(gca, 'FontSize', 12);
ylabel('Network number', 'FontSize', 20);
xlabel('Maximum $$t_{out}$$', 'FontSize', 20);

qsave = 0;
fname_str = 't_out_max_barh';
save_figure(h, 8, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

% display networks
networks_t10000 = find(max(t_out_all_net(:,:), [], 2)==tmax);
networks_t1000 = find(max(t_out_all_net(:,:), [], 2)>1000);

disp('tmax = 10000');
for i=1:numel(networks_t10000)
    disp( M_int_all_reduced{networks_t10000(i)} );
end

disp('tmax > 1000');
for i=1:numel(networks_t1000)
    disp( M_int_all_reduced{networks_t1000(i)} );
end

%% Mean t_out per topology

% selected networks to include
% default
network_sel_orig = network_all;
network_sel = 1:numel(network_all);
% custom
%network_sel_orig = network_set3_orig; %original index
%network_sel = network_set3; %reduced index

% get data
t_out_all_net = t_out_all(network_sel_orig, :, :);
hdata = mean(t_out_all_net(:,:), 2);

%{
figure;

histogram( hdata, 0:max(hdata)/20:max(hdata))
xlabel('Average t_{out}');
ylabel('Number of network topologies');
title('Highest period of different networks');
%}

% plot barh
h=figure;
hold on
[hdata_sorted, sort_idx] = sort(hdata, 'descend');
barh(hdata_sorted);
%plot([1 n_networks], [tmax tmax], 'r--');
plot([1 1], [0 numel(network_sel)+1], 'r-');
set(gca, 'YTick', 1:numel(network_sel), 'YTickLabels', string(network_sel(sort_idx)) );
set(gca, 'XScale', 'log');
xlim([10^(-0.1) 10^3]);
set(gca, 'FontSize', 12);
ylabel('Network number', 'FontSize', 20);
xlabel('Average $$t_{out}$$', 'FontSize', 20);

qsave = 0;
fname_str = 't_out_mean_barh';
save_figure(h, 8, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

% display networks
networks_t100 = find(mean(t_out_all_net(:,:), 2)>100);
networks_t40 = find(mean(t_out_all_net(:,:), 2)>40);

disp('max period > 100');
for i=1:numel(networks_t100)
    disp( M_int_all_reduced{networks_t100(i)} );
end

disp('max period > 40');
for i=1:numel(networks_t40)
    disp( M_int_all_reduced{networks_t40(i)} );
end

%% Look at some common topologies in more detail
% Histograms of mean t_out over different parameter sets
% Histograms of max t_out over different parameter sets
%{
%figure; 
edges = 0:tmax/20:tmax;
for i=1:numel(networks_t10000)
    %figure;
    %histogram( mean(t_out_all_net(networks_t10000(i), :, :), 3), edges );
    M_int = M_int_all_reduced{networks_t10000(i)};
    figure;
    histogram( max(t_out_all_net(networks_t10000(i), :, :), [], 3), edges );
    title(  sprintf('M_{int} = [%d %d; %d %d]', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2)) );
end
%}
%% -----------------Topologies capable of producing oscillations-----------
period_all_net = period_all(network_all, :, :);

[idx1, idx2] = find(period_all_net(:,:)~=Inf);
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
offset = signValues.*(width_change/2);
textPositions_cell = get(hText,{'Position'}); % cell array
textPositions = cell2mat(textPositions_cell); % numeric array
textPositions(:,1) = 1.3*textPositions(:,1) + offset; % add offset 

hText(1).Position = textPositions(1,:);
hText(2).Position = textPositions(2,:);

qsave = 0;
fname_str = 'oscillation_fraction_pie';
save_figure(h, 9, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Get list of topologies producing oscillations
topologies_osc = unique(idx1);
topologies_non_osc = setdiff(1:n_networks, unique(idx1));

%% Oscillations periods per topology
% selected networks to include
% default
%network_sel_orig = network_all;
%network_sel = 1:numel(network_all);
% subset
% network_sel_orig = network_set1_orig; %original index
% network_sel = network_set1; %reduced index
% custom
network_sel = [15 19 36 33 34 16 20 43]; %reduced index
network_sel_orig = network_all(network_sel); %original index

period_all_net = period_all(network_sel_orig, :, :);

h=figure;
hold on
% sort networks by number of different periods
%{
period_count = zeros(numel(network_sel), 1);
for i=1:numel(network_sel)
    period_count(i) = numel(unique(period_all_net(i,:)));
end
[~, idx_sort] = sort(period_count, 'descend');
for i=1:numel(network_sel)
    scatter( unique(period_all_net(idx_sort(i),:)), ...
        repmat(i, period_count(idx_sort(i)), 1), 35, 'kd', 'filled' );
end
set(gca, 'YTick', 1:numel(network_sel), 'YTickLabels', network_sel(idx_sort));

%}
% don't sort networks
for i=1:numel(network_sel)
    scatter( unique(period_all_net(i,:)), ...
        repmat(i, numel(unique(period_all_net(i,:))), 1), 35, 'kd', 'filled' );
end
set(gca, 'YTick', 1:numel(network_sel), 'YTickLabels', network_sel);

% set other properties
set(gca, 'XScale', 'log');
xlim([1 10^4]);
ylim([0 numel(network_sel)+1]);
colors = parula(4);
plot([2 2], [0 numel(network_sel)+1], '--', 'Color',  colors(1,:) );
plot([3 3], [0 numel(network_sel)+1], '--', 'Color', colors(2,:));
plot([4 4], [0 numel(network_sel)+1], '--', 'Color', colors(3,:));
set(gca, 'FontSize', 13);
ylabel('Network', 'FontSize', 20);
xlabel('Period', 'FontSize', 20);

% Same topologies as those with high tmax!
qsave = 0;
fname_str = 'oscillations_all_scatter_subset3_v2_reorganised';
save_figure(h, 10, 5, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Oscillations frequency per class (I, II, III)
frac_osc = zeros(3, 1); % Fraction of parameter sets which are capable of generating oscillations 
frac_osc_std = zeros(3, 1); % error of the mean across networks of same class
%frac_osc_2 = zeros(3, 1); % Net fraction of simulations which generate oscillations
for i=1:3
    period_temp = squeeze(period_all(network_sets{i}, :, :));
    n_network_temp = numel(network_sets{i});
    n_osc_psets = zeros(n_network_temp, 1);
    for j=1:n_network_temp
        [idx1, idx2] = find(squeeze(period_temp(j,:,:)~=Inf));
        n_osc_psets(j) = numel(unique(idx1));
    end
    
    frac_osc(i) = mean(n_osc_psets/n_pset); % count parameter sets
    frac_osc_std(i) = std(n_osc_psets/(n_pset) );
    %frac_osc_2(i) = numel(idx1)/(n_network_temp*n_pset*nsim); % count simulations
end

% Plot 
h = figure;
hold on
%bar([frac_osc([3 2 1]) 1-frac_osc([3 2 1])], 'stacked');
labels = {'Non-periodic\newline(n=9)', 'Oscillatory\newline(n=25)',...
    'Complex\newline(n=10)'};
set(gca, 'XTick', 1:3, 'XTickLabel', labels);
errorbar(frac_osc([3 2 1]), frac_osc_std([3 2 1]), 'LineWidth', 3);
%xlabel('Network class');
ylabel('Frequency');
%legend({'Oscillatory', 'Non-oscillatory'}, 'location', 'no');
set(gca, 'FontSize', 20);
set(gca, 'YTick', 0:0.2:1);
xlim([0.5 3.5]);
ylim([0 1]);

qsave = 1;
fname_str = 'oscillation_fraction_by_network_class_v4_errorbar_std';
save_figure(h, 10, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Categorize periodic trajectories
% Subdivide periodic trajectories into classes 
counts = cell(3,1);
cats = cell(3,1);
for i=1:3
    period_temp = squeeze(period_all(network_sets{i}, :, :));
    n_network_temp = numel(network_sets{i});
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

%% redistribute into fewer categories (by hand)
cats_new = {'2', '3', '4', '>4', 'Inf'};
counts_new = cell(2,1);
counts_new{1} = zeros(5, 1);
counts_new{1}(1) = counts{1}(1);
counts_new{1}(2) = counts{1}(2);
counts_new{1}(3) = counts{1}(3);
counts_new{1}(4) = sum(counts{1}(4:end-1));
counts_new{1}(5) = counts{1}(end);

counts_new{2} = zeros(5, 1);
counts_new{2}(1) = counts{2}(1);
counts_new{2}(3) = counts{2}(2);
counts_new{2}(5) = counts{2}(3);

% Plot counts of new categories
h1=figure;
%bar([counts{2}(1:end-1)/sum(counts{2}(1:end-1)); zeros(size(counts{2}(1:end-1)))], 'stacked')
%bar([zeros(size(counts{1}(1:end-1))); counts{1}(1:end-1)/sum(counts{1}(1:end-1))], 'stacked')
counts_normalized = cell(2,1);
counts_normalized{1} = counts_new{1}(1:end-1)/sum(counts_new{1}(1:end-1));
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
legend(cats_new(1:4), 'location', 'eo', 'FontSize', 20);

h2=figure;
counts_normalized{2} = counts_new{2}(1:end-1)/sum(counts_new{2}(1:end-1));
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
legend(cats_new(1:4), 'location', 'eo', 'FontSize', 20);

qsave = 1;
fname_str = 'oscillation_breakdown_complex_networks_pie_v2';
save_figure(h1, 10, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

qsave = 1;
fname_str = 'oscillation_breakdown_osc_networks_pie_v2';
save_figure(h2, 10, 8, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Oscillations frequency per network
period_all_net = period_all(network_all, :, :);

frac_osc = zeros(n_networks, 1); % Fraction of parameter sets which are capable of generating oscillations
frac_osc_2 = zeros(n_networks, 1); % Net fraction of simulations which generate oscillations
for i=1:size(period_all_net, 1)
    period_temp = squeeze(period_all_net(i,:,:));
    [idx1, idx2] = find(period_temp~=Inf);
    frac_osc(i) = numel(unique(idx1))/n_pset;
    frac_osc_2(i) = numel(idx1)/(n_pset*nsim);
end

%{
figure
hold on
plot(1:n_networks, frac_osc, 'bo');
plot(1:n_networks, frac_osc_2, 'ro');
xlabel('Topology number');
ylabel('Fraction periodic'); % fraction of parameter set which is capable of producing oscillations
%}
% rank-frequency plot
[frac_osc_sorted, sort_idx] = sort(frac_osc, 'descend');
h = figure;
%c = categorical(string(sort_idx)); string(sort_idx)
barh(frac_osc_sorted);
set(gca, 'YTick', 1:n_networks, 'YTickLabels', string(sort_idx));
set(gca, 'FontSize', 13);
ylabel('Network', 'FontSize', 20);
xlabel('Oscillation frequency', 'FontSize', 20);

disp('Most common networks:')
for i=1:5
    disp( M_int_all_reduced{ sort_idx(i) } );
end

disp('Non-oscillatory networks:')
idx = find(frac_osc==0);
for i=1:numel(idx)
    disp( M_int_all_reduced{ idx(i) } );
end

qsave = 0;
fname_str = 'oscillation_fraction_parameters_barh';
save_figure(h, 8, 10, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Oscillation ubiquity vs. interactions
osc_frac_vs_int = zeros(5);
count = zeros(5);
for i=1:n_networks
    M_int = M_int_all_reduced{i};
    n_act = sum(M_int(:)==1);
    n_rep = sum(M_int(:)==-1);
    count(n_act+1, n_rep+1) = count(n_act+1, n_rep+1)+1;
    osc_frac_vs_int(n_act+1, n_rep+1) = osc_frac_vs_int(n_act+1, n_rep+1) + frac_osc_2(i);
end
osc_frac_vs_int(count>0) = osc_frac_vs_int(count>0)./count(count>0);

n_int = 4;
h = figure;
%imagesc(count);
imagesc(0:n_int, 0:n_int, count);
[a, b] = meshgrid(0:n_int, 0:n_int);
s = string(count); s = reshape(s, 25, 1);
text(a(:), b(:), s, 'HorizontalAlignment', 'center',...
    'FontSize', 18, 'Color', 'w');
set(gca,'YDir', 'normal', 'XTick', 0:n_int, 'YTick', 0:n_int);
xlabel('Repressive interactions');
ylabel('Activating interactions');
set(gca, 'FontSize', 18);
xlim([-0.5 n_int+0.5]);
ylim([-0.5 n_int+0.5]);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
c = colorbar;
cmax = max(count(:));
colormap(viridis(cmax+1));
caxis([0 cmax+1]);
set(c, 'YTick', 1:6);
set(c, 'YTick', 0.5:1:cmax+0.5, 'YTickLabel', string(0:cmax) );
title('Number of circuits');

qsave = 0;
fname_str = 'network_count_vs_int_imagesc';
save_figure(h, 7, 6, fullfile(save_folder, fname_str),'.pdf', qsave);
%% 
h = figure;
%imagesc(count);
imagesc(0:n_int, 0:n_int, osc_frac_vs_int);
[a, b] = meshgrid(0:n_int, 0:n_int);
s = string(round(osc_frac_vs_int, 3)); s = reshape(s, 25, 1);
text(a(:), b(:), s, 'HorizontalAlignment', 'center',...
    'FontSize', 16, 'Color', 'w');
set(gca,'YDir', 'normal', 'XTick', 0:n_int, 'YTick', 0:n_int);
xlabel('Repressive interactions');
ylabel('Activating interactions');
set(gca, 'FontSize', 18);
xlim([-0.5 n_int+0.5]);
ylim([-0.5 n_int+0.5]);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
c = colorbar;
colormap(viridis);
caxis([0 1])
title('Oscillation prevalence');

qsave = 0;
fname_str = 'oscillation_frac_vs_int_imagesc';
save_figure(h, 7, 6, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Topologies showing complex dynamics
% Complexity defined as: (1) either period > 4, or (2) t_out > t_max (eq.
% not reached at t=t_max)
frac_complex_all = zeros(n_networks, 1);
for i=1:n_networks
    network = network_all(i);
    
    idx1 = find(period_all(network, :, :) > 4 & period_all(network, :, :) < Inf);
    idx2 = find(t_out_all(network, :,:)==tmax);
    idx12 = union(idx1, idx2);
    
    frac_complex_all(i) = numel(idx12)/(n_pset*nsim);
end

[frac_complex_all_sorted, sort_idx] = sort(frac_complex_all, 'descend');
idx_nonzero = find(frac_complex_all_sorted>0);
%%
h = figure;
barh(frac_complex_all_sorted(idx_nonzero) );
set(gca, 'YTick', 1:numel(idx_nonzero), 'YTickLabels', string(sort_idx(idx_nonzero)));
set(gca, 'FontSize', 13);
ylabel('Network', 'FontSize', 20);
xlabel('Complexity occurence', 'FontSize', 20);
set(gca, 'FontSize', 18);

qsave = 1;
fname_str = 'complexity_occurence_barh';
save_figure(h, 7, 6, fullfile(save_folder, fname_str),'.pdf', qsave);

%% Topologies capable of producing non-uniform equilibrium states
%frac_non_unif = zeros(n_networks,1);
%frac_non_unif = sum(non_uniform_all(network_all, :), 2); % fraction of simulations with non-uniform lattice
%frac_non_unif_2 = sum( sum(non_uniform_all(network_all, :, :), 3) > 0, 2)/n_pset;
frac_non_unif = sum( (non_uniform_all(network_all, :) > 1), 2)/(n_pset*nsim); % fraction of simulations with non-uniform lattice
frac_non_unif_2 = sum( sum(non_uniform_all(network_all, :, :), 3) > 1, 2)/n_pset; % fraction of parameter sets that give some non-uniform final state

figure;
hold on
plot(1:n_networks, frac_non_unif, 'bo');
plot(1:n_networks, frac_non_unif_2, 'rx');
xlabel('Topology number');
ylabel('Fraction non-uniform'); % fraction of parameter set which is capable of producing oscillations
ylim([0 1])
xlim([1 n_networks]);

%% Plot all circuit topologies
%{
for i=1:n_networks
    h = figure(7);
    plot_circuit(gz, a0, rcell, M_int_all_reduced{i}, [10 10], Coff, 10*ones(2), lambda(2) )

    qsave = 0;
    fname_str = sprintf('Network_%d', i);
    save_figure(h, 8, 4, fullfile(save_folder, 'all_networks', fname_str), '.pdf', qsave);
    close all;
end
%}
%}