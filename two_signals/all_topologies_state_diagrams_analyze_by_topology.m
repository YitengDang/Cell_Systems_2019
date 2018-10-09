%% Analyze statistics on all possible phases (up to equivalence)
% Do the analysis by grouping together results for each topology
clear variables
close all
clc
%%
% Settings
single_cell = 0;
sym = 0; % include symmetric partners?
save_folder_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\all_topologies';
save_folder_data = 'H:\My Documents\Multicellular automaton\data\two_signals\all_topologies';
%save_folder_fig = 'D:\Multicellularity\figures\two_signals\all_topologies'; % for figures
%save_folder_data = 'D:\Multicellularity\data\two_signals\all_topologies'; % for data
mastersave = 1; % switch off all save options

% Load data
load_path = 'H:\My Documents\Multicellular automaton\data\two_signals\all_topologies';
%load_path = 'D:\Multicellularity\data\two_signals\all_topologies';
labels = {'multi_cell', 'single_cell', 'multi_cell_all_incl', 'single_cell_all_incl'};
label = labels{single_cell+1 + 2*sym};
fname_str = sprintf('all_topologies_data_%s', label);

load(fullfile(load_path, fname_str));

% general 
n_data = numel(state_diagrams); % number of data points
n_genes = 2;
n_states = n_genes^2; % # states
n_int = n_genes^n_genes; % # interactions
%% Overview
% (1) number of distinct state diagrams
% (2) common state diagrams: (a) statistics, (b) plot of most common diagrams
% (3) Crude classification into (i) only ss, (ii) ss and oscillating, (iii)
% only osc
% (4) Distribution of number of steady states
% (5) Oscillation classification
% (6b) # steady states vs. # activation/repression interactions
% (6c) <osc. period> vs. # activation/repression interactions
% (6d) # oscillations vs. # steady states

%% (1) Statistics
% total number of distinct state diagrams
fprintf('Total number of diagrams: %d \n', n_data);
fprintf('Maximum number of distinct diagrams: %d \n', 2^(n_states^2));
unique_diagrams = {};
%unique_diagrams{1} = state_diagrams{1};
for i=1:n_data
    same = 0;
    for j=1:numel(unique_diagrams)
        if all(unique_diagrams{j}==state_diagrams{i})
            same = 1;
            break
        end
    end
    
    if ~same
        unique_diagrams{end+1} = state_diagrams{i};
    end
end
n_unique = numel(unique_diagrams);
fprintf('Found number of unique diagrams: %d \n', n_unique);
%% Get the network numbers of the list of M_int_all
% Make a list of M_int as in batch_sim_topologies_analyze
M = [0 1 -1]; % index to interaction
M_int_all_reduced = {}; 
done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
for k=1:3^4
    [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], k);
    gM = [i22 i21; i12 i11];
    M_int = [M(i11) M(i12); M(i21) M(i22)];
    if done(i11,i12,i21,i22)
    	continue
    elseif k==1
        continue
    else
        M_int_all_reduced{end+1} = M_int;
        done(i11,i12,i21,i22) = 1;
        done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
    end
end
n_networks = numel(M_int_all_reduced);

% Classify given M_int according to the given list
%{
M_int = [1 -1; 0 1];
idx = [];
for i=1:numel(M_int_all_reduced)
    if all(all( M_int_all_reduced{i} == M_int ))
        idx = i;
    end
end
%}
M_idx_all = zeros(n_data, 1);
M_idx_count = zeros(n_networks, 1); %number of counts per network
for i=1:n_data
    M_int = M_int_all{i};
    M_idx = 0;
    for j=1:numel(M_int_all_reduced)
        if all(all( M_int_all_reduced{j} == M_int ))
            M_idx = j;
            M_idx_count(M_idx) = M_idx_count(M_idx)+1;
            break
        end
    end
    M_idx_all(i) = M_idx;
end

%% (3) Crude classification
% -> Plot as pie chart
cat_all = zeros(n_networks, 2); % col 1: steady states present?, col 2: osc. present?
for i=1:n_data
    M_idx = M_idx_all(i);
    if M_idx==0
        continue
    end
    
    cat_all(M_idx, 1) = cat_all(M_idx, 1) || numel(steady_states{i})>0; 
    cat_all(M_idx, 2) = cat_all(M_idx, 2) || numel(cycles_all{i})>0;
    
    %{
    num_ss = numel(steady_states{i});
    num_cycles = numel(cycles_all{i});
    if num_ss>0
        cat_all(M_idx, 1) = 1;
    end
    if num_cycles>0
        cat_all(M_idx, 2) = 1;
    end
    %}
end

% categorization: (1) only steady states, (2) only cycles, (3) steady states & cycles,
cat = zeros(1, 3); 
cat_all_temp = cat_all*[1; 2];
for i=1:3
    cat(i) = sum( cat_all_temp == i);
end

% Check that all phases have been found
clc
fprintf('Total # phases found: %d \n', sum(cat));
fprintf('Total # phases: %d \n', n_data);

%figure;
%bar(1:3, cat);
%set(gca,'XTick', 1:3, 'XTickLabels', {'ss only','ss & cycles','cycles only'});
%xlabel('Steady states');
%ylabel('Count');
%set(gca, 'FontSize', 16);
%set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
%title(sprintf('Total: %d phases', n_data));

h=figure(1);
plotHandle = pie(cat);
idx_cat = find(cat~=0);
cmp = [1 0 0;  
    0 1 0;      
    0 0 1];
colormap(cmp(idx_cat, :));

hText = findobj(plotHandle, 'Type', 'text'); % text object handles
percentValues = get(hText, 'String'); % percent values
countValues = string(cat);
txt = {'ss: '; 'osc: '; 'ss and osc: '};
%combinedtxt = strcat(txt, percentValues);
combinedtxt = strcat(txt(idx_cat), percentValues);

combinedtxt2 = strcat( '(', countValues(idx_cat), ' networks)');
oldExtents_cell = get(hText, 'Extent'); 
if numel(idx_cat)>1
    oldExtents = cell2mat(oldExtents_cell); % numeric array
else
    oldExtents = oldExtents_cell;
end

%hText.String = {combinedtxt{3}, combinedtxt2{3}};
for i=1:numel(idx_cat)
    hText(i).String = {combinedtxt{i}, combinedtxt2{i}};
end
%hText(1).String = {combinedtxt{1}, combinedtxt2{1}};
%hText(2).String = {combinedtxt{2}, combinedtxt2{2}};
%hText(3).String = {combinedtxt{3}, combinedtxt2{3}};
set(plotHandle(2:2:end), 'FontSize', 14);

newExtents_cell = get(hText, 'Extent'); % cell array
if numel(idx_cat)>1
    newExtents = cell2mat(newExtents_cell); % numeric array 
else
    newExtents = newExtents_cell;
end

width_change = (newExtents(:,3)-oldExtents(:,3));
signValues = sign(oldExtents(:,1));
offset = signValues.*(width_change/2);
textPositions_cell = get(hText,{'Position'}); % cell array
textPositions = cell2mat(textPositions_cell); % numeric array
textPositions(:,1) = 1.3*textPositions(:,1) + offset; % add offset 

for i=1:numel(idx_cat)
    hText(i).Position = textPositions(i,:);
end
%hText(1).Position = textPositions(1,:);
%hText(2).Position = textPositions(2,:);
%hText(3).Position = textPositions(3,:);

% save figure
qsave = 1;
save_folder = save_folder_fig;
fname_str = sprintf('classification_types_by_topology_%s', label);
path_out = fullfile(save_folder, fname_str);
save_figure(h, 7, 5, path_out, '.pdf', qsave && mastersave)

%% (4) Number of stationary cell states
% N.B. for a single cell, this gives the actual # of fixed points, as
% transitions are unambiguous
n_steady_states_max = zeros(n_networks, 1);
n_steady_states_avg = zeros(n_networks, 1);

for i=1:n_data
    M_idx = M_idx_all(i);
    if M_idx==0
        continue
    end
    num = numel(steady_states{i});
    n_steady_states_avg( M_idx ) = n_steady_states_avg( M_idx ) + num;
    if num > n_steady_states_max(M_idx)
        n_steady_states_max(M_idx) = num;
    end
end
n_steady_states_avg = n_steady_states_avg./M_idx_count;

% histogram
close all
h2 = figure(2);
histogram( categorical(n_steady_states_max) );
xlabel('Max # stationary cell states');
ylabel('Count');
set(gca, 'FontSize', 16);
set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);
title(sprintf('Total: %d topologies', n_networks));

[n_steady_states_avg, sort_idx] = sort(n_steady_states_avg, 'descend');
h3 = figure(3);
barh(1:n_networks, n_steady_states_avg);
xlabel('Avg # stationary cell states');
ylabel('Network');
xlim([0 4])
set(gca, 'FontSize', 16, 'YTick', 1:n_networks, 'YTickLabels', sort_idx);
set(h3, 'Units', 'Inches', 'Position', [1 1 9 8]);
title(sprintf('Total: %d topologies', n_networks));

% save figures
qsave = 1;
save_folder = save_folder_fig;
fname_str = sprintf('num_steady_states_max_by_topology_%s', label);
path_out = fullfile(save_folder, fname_str);
save_figure(h2, 7, 6, path_out, '.pdf', qsave && mastersave)

qsave = 1;
save_folder = save_folder_fig;
fname_str = sprintf('num_steady_states_mean_by_topology_%s', label);
path_out = fullfile(save_folder, fname_str);
save_figure(h3, 8, 10, path_out, '.pdf', qsave && mastersave)

%% (6) heat map relations
% (a) distribution of phases according to # activators and # repressive
% (b) # steady states vs. activation/repression -> hypothesis: phases with more repressors
% have more oscillatory behaviour and fewer steady states.
% (c) fraction of oscillatory trajectories vs. activating & repressive interactions
% (d) # oscillation period vs. activation/repression -> hypothesis: phases with more repressors
% have more oscillatory behaviour and fewer steady states.

% ------Calculations (process data)---------
count_int_type = zeros(n_int+1, n_int+1); % number of topologies with given interactions
for i=1:n_networks
    M_int = M_int_all_reduced{i};
    n_act = sum(M_int(:)==1); % number of activating interactions
    n_rep = sum(M_int(:)==-1); % number of repressing interactions
    count_int_type(n_act+1, n_rep+1) = count_int_type(n_act+1, n_rep+1) + 1;
end

ss_int_type = zeros(n_int+1, n_int+1, n_networks); % # steady states 
osc_q_int_type = zeros(n_int+1, n_int+1, n_networks); % oscillatory?
osc_period_int_type = zeros(n_int+1, n_int+1, n_networks); % average oscillation period

%C = categorical(cycle_structures);
%osc_vs_ss = zeros(n_states+1, numel(unique(C))); % #osc. vs. #steady states

for i=1:n_data
    M_int = M_int_all{i};
    M_idx = M_idx_all(i);
    if M_idx==0
        continue
    end
    n_act = sum(M_int(:)==1); % number of activating interactions
    n_rep = sum(M_int(:)==-1); % number of repressing interactions
    n_ss = numel(steady_states{i}); % number of steady states
        
    ss_int_type(n_act+1, n_rep+1, M_idx) = ss_int_type(n_act+1, n_rep+1, M_idx) + n_ss;
    %
    if ~isempty(cycles_all{i})
        osc_q_int_type(n_act+1, n_rep+1, M_idx) = ...
            osc_q_int_type(n_act+1, n_rep+1, M_idx) + 1;
        av_osc = mean( cellfun(@numel, cycles_all{i})-1 );
        osc_period_int_type(n_act+1, n_rep+1, M_idx) = ...
            osc_period_int_type(n_act+1, n_rep+1, M_idx) + av_osc;
    end
    %}
    %idx_C = find(C(i) == categories);
    %idx_C = find(C(i) == cat_sorted);
    %osc_vs_ss(n_ss+1, idx_C) = osc_vs_ss(n_ss+1, idx_C) + 1;
    %scatter(n_ss, n_int, 'bo');
end

% calculate average # steady states
idx = (count_int_type~=0);

% take mean over all networks
for i=1:n_networks
    ss_int_type(:,:,i) = ss_int_type(:,:,i)/M_idx_count(i);
    
    osc_q_temp = sum(sum(osc_q_int_type(:,:,i))); % number of oscillatory phases for given network
    osc_q_int_type(:,:,i) = osc_q_int_type(:,:,i)/M_idx_count(i);
    if osc_q_temp~=0
        osc_period_int_type(:,:,i) = osc_period_int_type(:,:,i)/osc_q_temp;        
    end
end
ss_data = sum( ss_int_type, 3 );
ss_data(idx) = ss_data(idx)./count_int_type(idx);

osc_q_data = sum( osc_q_int_type, 3 );
osc_q_data(idx) = osc_q_data(idx)./count_int_type(idx);

osc_period_data = sum( osc_period_int_type, 3 );
osc_period_data(idx) = osc_period_data(idx)./count_int_type(idx);
%% -----------Plots--------------
% Plot (a)
h = figure(6);
hold on
imagesc(0:n_int, 0:n_int, count_int_type);
[a, b] = meshgrid(0:n_int, 0:n_int);
s = string(count_int_type); s = reshape(s, 25, 1);
text(a(:), b(:), s, 'HorizontalAlignment', 'center',...
    'FontSize', 14, 'Color', 'w');
set(gca,'YDir', 'normal', 'XTick', 0:n_int, 'YTick', 0:n_int);
xlabel('Repressive interactions');
ylabel('Activating interactions');
set(gca, 'FontSize', 18);
xlim([-0.5 n_int+0.5]);
ylim([-0.5 n_int+0.5]);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
c = colorbar;
cm_viridis = viridis;
colormap(cm_viridis);
ylabel(c, 'Count');
title(sprintf('Number of phases (total: %d)', n_data));

% save figure
qsave = 1;
save_folder = save_folder_fig;
fname_str = sprintf('act_rep_heatmap_topology_count_%s', label);
path_out = fullfile(save_folder, fname_str);
save_figure(h, 7, 6, path_out, '.pdf', qsave && mastersave)
%% (b) # steady states vs. activation/repression
h = figure(7);
hold on
imagesc(0:n_int, 0:n_int, ss_data);
[a, b] = meshgrid(0:n_int, 0:n_int);
s = string(round(ss_data, 3)); s = reshape(s, 25, 1);
text(a(:), b(:), s, 'HorizontalAlignment', 'center',...
    'FontSize', 14, 'Color', 'w');
set(gca,'YDir', 'normal', 'XTick', 0:n_int, 'YTick', 0:n_int);
xlabel('Repressive interactions');
ylabel('Activating interactions');
set(gca, 'FontSize', 18);
xlim([-0.5 n_int+0.5]);
ylim([-0.5 n_int+0.5]);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
c = colorbar;
cm_viridis = viridis;
colormap(cm_viridis);
caxis([0 n_states])
title('Mean number of stationary cell states');

% save figure
qsave = 1;
save_folder = save_folder_fig;
fname_str = sprintf('act_rep_heatmap_num_ss_%s', label);
path_out = fullfile(save_folder, fname_str);
save_figure(h, 7, 6, path_out, '.pdf', qsave && mastersave)

%% (c) fraction of oscillatory trajectories vs. activating & repressive interactions
h = figure;
hold on
imagesc(0:n_int, 0:n_int, osc_q_data);
[a, b] = meshgrid(0:n_int, 0:n_int);
s = string(round(osc_q_data, 3)); s = reshape(s, numel(s), 1);
text(a(:), b(:), s, 'HorizontalAlignment', 'center',...
    'FontSize', 14, 'Color', 'w');
set(gca,'YDir', 'normal', 'XTick', 0:n_int, 'YTick', 0:n_int);
xlabel('Repressive interactions');
ylabel('Activating interactions');
set(gca, 'FontSize', 18);
xlim([-0.5 n_int+0.5]);
ylim([-0.5 n_int+0.5]);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
c = colorbar;
cm_viridis = viridis;
colormap(cm_viridis);
ylabel(c, 'Fraction');
caxis([0 1]);
title('Diagrams with cycles');

% save figure
qsave = 1;
save_folder = save_folder_fig;
fname_str = sprintf('act_rep_heatmap_frac_osc_%s', label);
path_out = fullfile(save_folder, fname_str);
save_figure(h, 7, 6, path_out, '.pdf', qsave && mastersave)

%% (d) oscillation period vs. activating/repressing interactions
h = figure(9);
hold on
imagesc(0:n_int, 0:n_int, osc_period_data);
[a, b] = meshgrid(0:n_int, 0:n_int);
s = string(round(osc_period_data, 3)); s = reshape(s, numel(s), 1);
text(a(:), b(:), s, 'HorizontalAlignment', 'center',...
    'FontSize', 14, 'Color', 'w');
set(gca,'YDir', 'normal', 'XTick', 0:n_int, 'YTick', 0:n_int);
xlabel('Repressive interactions');
ylabel('Activating interactions');
set(gca, 'FontSize', 18);
xlim([-0.5 n_int+0.5]);
ylim([-0.5 n_int+0.5]);
set(h, 'Units', 'Inches', 'Position', [1 1 9 8]);
c = colorbar;
cm_viridis = viridis;
colormap(cm_viridis);
ylabel(c, 'Fraction');
caxis([0 4]);
title('Average cycle length');

% save figure
qsave = 1;
save_folder = save_folder_fig;
fname_str = sprintf('act_rep_heatmap_osc_period_%s', label);
path_out = fullfile(save_folder, fname_str);
save_figure(h, 7, 6, path_out, '.pdf', qsave && mastersave)

%% Functions
function h = draw_state_diagram(A, fig_num, count)
    % A: state diagram (graph adjacency matrix)
    % fig_num: number of figure to plot
    % count: number of times this diagram appears among all phases
    % h: returns figure handle
    if nargin<2
        fig_num = 1;
        count = 0;
    end
    h = figure(fig_num);
    hold on
    s = [0 1 0 1];
    t = [0 0 1 1];
    %A = ones(4);
    Gs = digraph(A);
    nLabels = {};
    g=plot(Gs, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
        'LineWidth', 3, 'EdgeColor', 'k',...
        'Marker', 'o', 'MarkerSize', 100, 'NodeColor', [0.2 0.2 0.2], 'NodeLabel', nLabels);
    % Make edges dashed if state has two outgoing edges
    for i=1:4
        if sum(A(i,:))==2
            idx = find(A(i,:));
            highlight(g, i, idx, 'LineStyle', '--', 'LineWidth', 2);
        elseif sum(A(i,:))==4
            highlight(g, i, 1:4, 'LineStyle', ':', 'LineWidth', 2);
        end
    end
    text(s-0.11,t+0.019,{'(0,0)','(1,0)','(0,1)','(1,1)'}, 'Color', 'w', 'FontSize', 32)
    
    text(0.5, 1.3, sprintf('Number of diagrams = %d', count), ...
        'FontSize', 16, 'HorizontalAlignment', 'center')
    ax = gca;
    axis([-0.4 1.4 -0.4 1.4]);
    ax.Visible = 'off';
    h.Color = [1 1 1];
    set(ax, 'Units', 'Inches', 'Position', [0 0 7 6]);
    set(h, 'Units', 'Inches', 'Position', [0.2 0.2 7 6]);
end