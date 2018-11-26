% Tests whether a travelling wave can propagate by performing explicit
% simulations
clear all
close all
set(0,'defaulttextinterpreter', 'latex')
%% Parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];

Con = [18 16];
Coff = [1 1];
M_int = [0 1; -1 1];
K = [0 9; 11 6];

hill = Inf;
noise = 0;

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Specify wave type
num_waves = 1;
wave_type = 3;
fname_str_all = {'trav_wave_single_vertical',...
    'trav_wave_single_horizontal_inward_bend',...
    'trav_wave_single_horizontal_outward_bend'};
fname_str = fname_str_all{wave_type};

%% Pre-prosessing
% calculate fN
idx_loop = gz + round(gz/2); % pick cell not at corner of grid
dist_vec = a0*dist(idx_loop,:);
r = dist_vec(dist_vec>0);
fN = zeros(2,1);
fN(1) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
fN(2) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));

% save folder
%save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_stability';
%fname_str_default = strrep(sprintf('N%d_a0_%.1f_rcell_%.1f_lambda12_%.1f_M_int_%d_%d_%d_%d_Con_%d_%d',...
%    N, a0, rcell, lambda(2), M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2),...
%    Con(1), Con(2)), '.', 'p');

% Load initial conditions
%{
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_snapshots';

fname = fullfile(load_folder, fname_str);
cells_load = cell(2,1);
cells_load{1} = xlsread(fname, 'Sheet1');
cells_load{2} = xlsread(fname, 'Sheet2');
if all(size(cells_load{1})==[N 2])
    cells_in = cells_load{1};
elseif all(size(cells_load{1})==[gz gz]) && all(size(cells_load{2})==[gz gz])
    cells_in(:, 1) = reshape(cells_load{1}, N, 1);
    cells_in(:, 2) = reshape(cells_load{2}, N, 1);
else
    disp('Wrong input format');
end

% get cell state population of initial state
cells_idx = cells_in*[1; 2];
n_in = histcounts(cells_idx, -0.5:3.5);
%}

%%
% Load analyzed data
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
subfolder = 'run2';
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_%s',...
    num_waves, wave_type, subfolder);
load( fullfile(load_folder, fname_str) );
%}

% display found waves
%
[x_found, y_found] = find(wave_possible);
num_psets = numel(x_found);
t=table(P(x_found, :), y_found, 'VariableNames', {'Wave_type', 'Network'});
t2=table(x_found, y_found, 'VariableNames', {'wave_idx', 'Network'});

disp('Cell states F, M, B, E');
disp(t);
disp(t2);

% Get a list of all interaction matrices for the found networks
%[M_int_found] = get_found_M_int(y_found);

% Load new data from reanalyzed wave_type=2 and wave_type=3
% old code
%{
num_psets = numel(x_found);
frac_correct = zeros(num_psets, 1);
frac_correct_strict = zeros(num_psets, 1);
num_predicted_all = zeros(num_psets, 1);
%}
% output variables
results_all = cell(num_psets, 1); % aggregate results 
% out_all: for each parameter set, determine whether it is a (1) true
% negative, (2) false negative, (3) false positive, (4) true positive
out_all_count = zeros(num_psets, 4);

for idx_loop=1:num_psets
    %%
    %idx = 1;
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    states_perm = P(wave_idx, :);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    
    % Load rerun data
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2_rerun';
    if wave_type==1
        fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d', ...
            num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4) );
    elseif wave_type==2 || wave_type==3
        fname_str = sprintf('trav_wave_conditions_rerun_all_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d',...
            num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    end
    fname = fullfile(load_folder, fname_str);
    load( fname, 'trav_wave_conds_met_both_types' );
    %%
    % Load data from simulations 
    labels = {'single_vertical', 'single_horizontal_inward_bend', 'single_horizontal_outward_bend'};    
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2_stability_sim';
    fname_str = sprintf(...
        'stability_sim_from_pred_all_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d',...
        num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4) );
    load( fullfile(load_folder, fname_str), 'trav_wave_all', 'unmodified_all' );
    
    trav_wave_unmodified_all = trav_wave_all & unmodified_all;
    %%
    % store results
    results_all{idx_loop}(:, 1) = trav_wave_conds_met_both_types; % column 1: theory
    results_all{idx_loop}(:, 2) = trav_wave_unmodified_all; % column 2: simulations
    out_all = results_all{idx_loop}*[2; 1]; 
    for i=1:4
        out_all_count(idx_loop, i) = sum(out_all == i-1);
    end
    
    % Old code
    %{
    % Get list of which parameter sets should be included
    n_data_raw = size(Con_wave, 1);
    pset_idx = zeros(n_data_raw, 1);
    for i=1:n_data_raw
        pset_idx(i) = ~isempty(find(Con_wave(i,:)==Con_wave_new, 1));
    end
    
    % for these included parameter sets, determine how many are correctly
    % predicted
    n_predicted = size(Con_wave_new, 1); % number of trav. wave parameter sets predicted, both conditions satisfied
    n_actual = sum(trav_wave_all & pset_idx);
    n_actual_strict = sum(trav_wave_all_strict & pset_idx);
    
    frac_correct(idx) = n_actual/n_predicted;
    frac_correct_strict(idx) = n_actual_strict/n_predicted;
    num_predicted_all(idx) = n_predicted;
    %}
end
%% Plot overall density of travelling waves (from simulations)
Qvals_all = out_all_count(:, 4)/size(results_all{1}, 1);

h = figure;
hold on
bar(Qvals_all);

% labels below
%labels = {'15, [3 4 2 1]', '19, [4 3 1 2]', '33, [4 2 1 3]', '33, [3 4 2 1]',...
%    '34, [4 2 1 3]', '36, [2 4 3 1]'};
labels = {'15', '19', '33(a)', '33(b)',...
    '34', '36'};
set(gca, 'XTick', 1:num_psets, 'XTickLabels', labels, 'XTickLabelRotation', 45);

% other settings
%xlabel('Wave type');
xlabel('Network');
ylabel('Robustness');
set(gca, 'FontSize', 20);
set(h, 'Units', 'inches', 'position', [1 1 10 6]);
ylim([0 0.02]);
%}

qsave = 1;    
pred_label = 'run2';
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\all_networks_test_analytical_in_sims';
fname_str_save = sprintf('wave_num_%d_wave_type_%d_%s_Q_vals', num_waves, wave_type, pred_label);
fname = fullfile(save_folder, fname_str_save);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%% normalize by the number of parameters for each circuit
num_params = [5 5 6 6 6 6];
Qvals_all_norm = (Qvals_all').^(1./num_params);

h = figure;
hold on
bar(Qvals_all_norm);

% labels below
%labels = {'15, [3 4 2 1]', '19, [4 3 1 2]', '33, [4 2 1 3]', '33, [3 4 2 1]',...
%    '34, [4 2 1 3]', '36, [2 4 3 1]'};
labels = {'15', '19', '33(a)', '33(b)',...
    '34', '36'};
set(gca, 'XTick', 1:num_psets, 'XTickLabels', labels, 'XTickLabelRotation', 45);

% other settings
%xlabel('Wave type');
xlabel('Network');
ylabel('Robustness (normalized)');
set(gca, 'FontSize', 20);
set(h, 'Units', 'inches', 'position', [1 1 10 6]);
ylim([0 1]);
%}

qsave = 1;    
pred_label = 'run2';
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\all_networks_test_analytical_in_sims';
fname_str_save = sprintf('wave_num_%d_wave_type_%d_%s_Q_vals_normalized_2', num_waves, wave_type, pred_label);
fname = fullfile(save_folder, fname_str_save);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Calculate precision, recall, F1 score
precision_all = out_all_count(:, 4)./(out_all_count(:, 3)+out_all_count(:, 4));
recall_all = out_all_count(:, 4)./(out_all_count(:, 2)+out_all_count(:, 4));
F1_score_all = 2*precision_all.*recall_all./(precision_all + recall_all);
disp('precision_all');
disp(precision_all);
disp('recall_all');
disp(recall_all);
disp('F1_score_all');
disp(F1_score_all);

%% Plot P, R, F1 score together
h = figure;
hold on
bar_data = [precision_all recall_all F1_score_all];
%bar_data = [precision_all recall_all];
bar(bar_data)

% labels below
%labels = {'15, [3 4 2 1]', '19, [4 3 1 2]', '33, [4 2 1 3]', '33, [3 4 2 1]',...
%    '34, [4 2 1 3]', '36, [2 4 3 1]'};
labels = {'15', '19', '33(a)', '33(b)',...
    '34', '36'};
set(gca, 'XTick', 1:num_psets, 'XTickLabels', labels, 'XTickLabelRotation', 45);

% labels above: number of true positives
labels = sprintfc('TP=%d', out_all_count(:, 4));
xt = get(gca, 'XTick');
text(xt, max(bar_data, [], 2), labels, 'FontSize', 20,...
    'HorizontalAlignment', 'center', 'VerticalAlignment','bottom')

% other settings
%xlabel('Wave type');
xlabel('Network');
ylabel('Value');
set(gca, 'FontSize', 20);
set(h, 'Units', 'inches', 'position', [1 1 10 6]);
ylim([0 1]);
legend({'Precision', 'Recall', 'F1 score'}, 'Location', 'ne');
%}

qsave = 1;    
pred_label = 'run2_rerun';
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\all_networks_test_analytical_in_sims';
fname_str_save = sprintf('wave_num_%d_wave_type_%d_%s_P_R_F1_together', num_waves, wave_type, pred_label);
%fname_str_save = sprintf('wave_num_%d_wave_type_%d_%s_P_R_together', num_waves, wave_type, pred_label);
fname = fullfile(save_folder, fname_str_save);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%% Plot precision
h = figure;
hold on
bar(precision_all);

% labels below
%labels = {'15, [3 4 2 1]', '19, [4 3 1 2]', '33, [4 2 1 3]', '33, [3 4 2 1]',...
%    '34, [4 2 1 3]', '36, [2 4 3 1]'};
labels = {'15', '19', '33(a)', '33(b)',...
    '34', '36'};
set(gca, 'XTick', 1:num_psets, 'XTickLabels', labels, 'XTickLabelRotation', 45);

% labels above
%{
labels = sprintfc('n=%d', num_predicted_all);
xt = get(gca, 'XTick');
text(xt, frac_correct, labels, 'FontSize', 20,...
    'HorizontalAlignment', 'center', 'VerticalAlignment','bottom')
%}

% other settings
%xlabel('Wave type');
xlabel('Network');
ylabel('Precision');
set(gca, 'FontSize', 20);
set(h, 'Units', 'inches', 'position', [1 1 10 6]);
ylim([0 1]);
%}

qsave = 1;    
pred_label = 'run2_rerun';
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\all_networks_test_analytical_in_sims';
fname_str_save = sprintf('wave_num_%d_wave_type_%d_%s_precision', num_waves, wave_type, pred_label);
fname = fullfile(save_folder, fname_str_save);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot recall
h = figure;
hold on
bar(recall_all);

% labels below
%labels = {'15, [3 4 2 1]', '19, [4 3 1 2]', '33, [4 2 1 3]', '33, [3 4 2 1]',...
%    '34, [4 2 1 3]', '36, [2 4 3 1]'};
labels = {'15', '19', '33(a)', '33(b)',...
    '34', '36'};
set(gca, 'XTick', 1:num_psets, 'XTickLabels', labels, 'XTickLabelRotation', 45);

% other settings
%xlabel('Wave type');
xlabel('Network');
ylabel('Recall');
set(gca, 'FontSize', 20);
set(h, 'Units', 'inches', 'position', [1 1 10 6]);
ylim([0 1]);
%}

qsave = 1;    
pred_label = 'run2_rerun';
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\all_networks_test_analytical_in_sims';
fname_str_save = sprintf('wave_num_%d_wave_type_%d_%s_recall', num_waves, wave_type, pred_label);
fname = fullfile(save_folder, fname_str_save);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot F1 score
h = figure;
hold on
bar(F1_score_all);

% labels below
%labels = {'15, [3 4 2 1]', '19, [4 3 1 2]', '33, [4 2 1 3]', '33, [3 4 2 1]',...
%    '34, [4 2 1 3]', '36, [2 4 3 1]'};
labels = {'15', '19', '33(a)', '33(b)',...
    '34', '36'};
set(gca, 'XTick', 1:num_psets, 'XTickLabels', labels, 'XTickLabelRotation', 45);

% other settings
%xlabel('Wave type');
xlabel('Network');
ylabel('F1 score');
set(gca, 'FontSize', 20);
set(h, 'Units', 'inches', 'position', [1 1 10 6]);
ylim([0 1]);
%}

qsave = 1;    
pred_label = 'run2_rerun';
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\all_networks_test_analytical_in_sims';
fname_str_save = sprintf('wave_num_%d_wave_type_%d_%s_F1_score', num_waves, wave_type, pred_label);
fname = fullfile(save_folder, fname_str_save);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot correct and wrong data points together on phase maps
% Load analyzed data
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
subfolder = 'run2';
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_%s',...
    num_waves, wave_type, subfolder);
load( fullfile(folder, fname_str) );

% get waveforms and networks of possible waves
[x_found, y_found] = find(wave_possible);

% Get a list of all interaction matrices for the found networks
%[M_int_found] = get_found_M_int(y_found);

for idx_loop=1:numel(x_found)
    %idx_loop = 1;
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    states_perm = P(wave_idx, :);
    
    % -----Load data... same as previous-----------------------------------
    % Load original data 
    labels = {'single_vertical', 'single_horizontal_inward_bend', 'single_horizontal_outward_bend'};    
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
    fname_str = sprintf(...
        'stability_sim_from_pred_trav_wave_%s_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d',...
        labels{wave_type}, num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4) );
    load( fullfile(load_folder, fname_str) );
    %trav_wave_all = trav_wave_conds_met(network, :);
    
    % Load rerun data
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2_rerun';
    fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d', ...
        num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4) );
    load( fullfile(load_folder, fname_str) );
    
    % Get list of which parameter sets should be included
    n_data_raw = size(Con_wave, 1);
    pset_idx = zeros(n_data_raw, 1);
    for i=1:n_data_raw
        pset_idx(i) = ~isempty(find(Con_wave(i,:)==Con_wave_new, 1));
    end
    %----------------------------------------------------------------------
    % Filter data that has all travelling wave
    % conditions satisfied (e.g. for a type 2 wave, the rerun data gives
    % whether the type 1 conditions are also satisfied).
    
    idx = find(pset_idx);
    Con_wave_filtered = Con_wave(idx, :);
    K_wave_filtered = K_wave(idx, :, :);
    trav_wave_filtered = trav_wave_all(idx);
    
    % Then, distinguish between true positives (theory & simulations ->
    % trav. wave) and false positives (theory-> TW, simulations-> not TW)
    idx_true = find(trav_wave_filtered);
    idx_false = find(~trav_wave_filtered);
    Con_wave_true = Con_wave(idx_true, :);
    Con_wave_false = Con_wave(idx_false, :);
    K_wave_true = K_wave(idx_true, :, :);
    K_wave_false = K_wave(idx_false, :, :);
    %}
    %--------------------------------------------------------------------------
    %% Plot phase diagram
    %Con_wave = squeeze(Con_all(network, squeeze(trav_wave_conds_met(network, :)), :));
    %K_wave = squeeze(K_all(network, squeeze(trav_wave_conds_met(network, :)), :, :));
    for idx_i=1:2
        for idx_j=1:2
            %idx_i = 1;
            %idx_j = 1;

            % check that the interaction is not absent
            if K_wave(1, idx_i, idx_j)==0
                continue
            end

            h = plot_phase_diagram_local(a0, rcell, lambda, dist, idx_i, idx_j);
            %h=figure;
            hold on
            color_true = [0.1 0.1 0.1];
            p1=scatter(K_wave_true(:,idx_i,idx_j), Con_wave_true(:,idx_j), 'o',...
                'MarkerEdgeColor', color_true, 'MarkerFaceColor', color_true, 'MarkerFaceAlpha', 0.5);
            color_false = [0.9 0.6 0.9];
            p2=scatter(K_wave_false(:,idx_i,idx_j), Con_wave_false(:,idx_j), 'o',...
                'MarkerEdgeColor', color_false, 'MarkerFaceColor', color_false, 'MarkerFaceAlpha', 0.5);
            %scatter([500 510 520 530], [100 90 80 110], 'x', 'MarkerEdgeColor', [0.9 0.6 0.9]);
            legend([p1, p2], {'True positive', 'False positive'}, 'Location', 'se', 'FontSize', 16);
            
            % save figures
            qsave = 1;
            save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\all_networks_analytical_phases\phase_diagram';
            fname_str = sprintf('Trav_wave_conditions_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d_interaction_%d_%d_scatter_sim_tested',...
                num_waves, wave_type, network, states_perm(1), states_perm(2),...
                states_perm(3), states_perm(4), idx_i, idx_j);
            fname = fullfile(save_folder, 'sim_tested', fname_str);
            save_figure(h, 10, 8, fname, '.pdf', qsave);
        end
    end
end

%% Functions
function h = plot_phase_diagram_local(a0, rcell, lambda, dist, idx_gene, idx_mol)
    % idx_gene: index of gene under control
    % idx_mol: index of sensed molecule / molecule affecting the gene (1 <= i <= L)
    
    % calculate fN
    Rcell = rcell*a0;
    fN = zeros(2,1);
    %[dist, ~] = init_dist_hex(gz, gz);
    dist_vec = a0*dist(1,:);
    r = dist_vec(dist_vec>0); % exclude self influence
    for i=1:2 % calculate signaling strength
        fN(i) = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(i)).*(lambda(i)./r)) ); 
    end

    % Parameter range map
    num_points = 1000;
    Con_max = 1000;
    K_max = 1000;
    Con_vec = linspace(1, Con_max, num_points);
    K_vec = linspace(1, K_max, num_points);
    [K, Con] = meshgrid(K_vec, Con_vec);

    % Make 4 limiting regions as boolean matrices
    R1 = (1+fN(idx_mol) - K) > 0; % Everything ON
    R2 = ((1+fN(idx_mol))*Con - K ) < 0; % Everything OFF
    R3 = ((Con + fN(idx_mol) - K) > 0 & (1+fN(idx_mol) - K) < 0); % ON remains ON & not all ON
    R4 = ((1+ fN(idx_mol)*Con - K) < 0 & ((1+fN(idx_mol))*Con - K ) > 0) ; % OFF remains OFF & not all OFF
    %R3 = (Con > K & (K-1)./Con > fN); % autonomous cells for Son > K
    %R4 = (Con <= K & K - Con < fN & (K-1)./Con > fN); % autonomous cells for Son < K

    out = R1 + 2*R2 + 3*R3 + 4*R4; % only regions 3 and 4 can overlap
    if ~isempty(find(unique(out)==0, 1))
        map_idx = 5; % activation-deactivation
        out(out==0) = 5; 
        %phase = 'non-A';
        phase = 'U';
    elseif ~isempty(find(unique(out)==7, 1)) 
        map_idx = 6; % autonomy
        out(out==7) = 5; % ON remains ON & OFF remains OFF
        phase = 'A01';
    end

    % plot figure
    h = figure;
    hold on
    him = imagesc(K_vec, Con_vec, out);

    %set(him, 'AlphaData', out > 0); % invisible if not from any region
    % R1 -> black
    % R2 -> white
    % R3 -> green
    % R4 -> red
    % activation-deactivation -> magenta
    % autonomy -> gray
    map = [0, 0, 0
        1, 1, 1
        0, 1, 0
        1, 0, 0
        1, 1, 0
        0.5, 0.5, 0.5];
    tmp = map([1:4 map_idx], :);
    colormap(tmp);
    c=colorbar;
    set(c, 'YTick', 1+2/5+4/5*(0:4));

    %phase_labels = {'all>K','all<K','A1','A0',phase};
    phase_labels = {'P1','P0','A1','A0',phase};
    set(c, 'TickLabels', phase_labels);

    % adjust the graph
    set(gca,'ydir', 'normal', 'FontSize', 24)
    xlabel('K', 'FontSize', 24)
    ylabel('$$C_{ON}$$', 'FontSize', 24)
    ylim([1 Con_max])
    xlim([1 K_max])
    %title(sprintf('$$f_N = %.3f, a_0 = %.2f$$', fN, a0), 'FontSize', 30)
    title(sprintf('Interaction $$%d \\leftarrow %d$$', idx_gene, idx_mol))
    %title(sprintf('Molecule %d', idx_mol));
end
