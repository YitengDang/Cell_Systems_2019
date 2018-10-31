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
wave_type = 1;
fname_str_all = {'trav_wave_single_vertical',...
    'trav_wave_single_horizontal_inward_bend',...
    'trav_wave_single_horizontal_outward_bend'};
fname_str = fname_str_all{wave_type};

%% Pre-prosessing
% calculate fN
idx = gz + round(gz/2); % pick cell not at corner of grid
dist_vec = a0*dist(idx,:);
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

%% Load data
% v1 data
%{
load_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\data';
%fname_str_data = sprintf('%s_stability_sim_%s_K12_K21_K22_range', fname_str, fname_str_default);  
fname_str_data = sprintf('%s_stability_sim_K12_K21_K22_range', fname_str);  
fname = fullfile(load_folder, fname_str_data);
load(fname, 'K_12_all', 'K_21_all', 'K_22_all', 'trav_wave_all', 'unmodified_all');
%}
% v2 data
load_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\data';
fname_str_temp = {'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_15_states_3_4_2_1',...
'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_19_states_4_3_1_2',...
'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_33_states_3_4_2_1',...
'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_33_states_4_2_1_3',...
'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_34_states_4_2_1_3',...
'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_36_states_2_4_3_1'};
fname_str_all = cellfun(@(x) sprintf(x, num_waves, wave_type), fname_str_temp, 'UniformOutput', false);

num_files = numel(fname_str_all);
frac_correct = zeros(num_files, 1);
frac_correct_strict = zeros(num_files, 1);
num_sets_all = zeros(num_files, 1);

for idx=1:num_files
    fname = fullfile(load_folder, fname_str_all{idx});
    load(fname);
    
    frac_correct(idx) = sum(trav_wave_all)/numel(trav_wave_all);
    frac_correct_strict(idx) = sum(trav_wave_all_strict)/numel(trav_wave_all_strict);
    num_sets_all(idx) = numel(trav_wave_all);
end
%}
%{
trav_wave_temp1 = trav_wave_all;
trav_wave_temp2 = trav_wave_all_strict;
trav_wave_all = trav_wave_temp2;
trav_wave_all_strict = trav_wave_temp1;
save(fname, 'M_int', 'fname_str_pred', 'network', 'states_perm',...
    'trav_wave_all_strict', 'trav_wave_all', 'unmodified_all', 'wave_idx');
%}
%% Plot fraction of correct predictors
h = figure;
hold on
bar(frac_correct);
%bar([frac_correct, frac_correct_strict]);
%xlim([K_12_all(1)-1 K_12_all(end)+1]);
%ylim([K_21_all(1)-1 K_21_all(end)+1]);
labels = {'15, [3 4 2 1]','19, [4 3 1 2]','33, [4 2 1 3]', '33, [3 4 2 1]',...
    '34, [4 2 1 3]','36, [2 4 3 1]'};
set(gca, 'XTick', 1:num_files, 'XTickLabels', labels, 'XTickLabelRotation', 45);
xlabel('Wave type');
ylabel('Fraction predicted correctly');
set(gca, 'FontSize', 20);
set(h, 'Units', 'inches', 'position', [1 1 10 6]);
ylim([0 1]);
%}

qsave = 0;    
pred_label = 'run2';
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\all_networks_test_analytical_in_sims';
fname_str_save = sprintf('frac_correctly_predicted_by_network_wave_num_%d_wave_type_%d_%s', num_waves, wave_type, pred_label);
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
    
    % load NNMFA data
    %folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
    %subfolder = 'run3_vary_a0_lambda12';
    fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
           num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    %disp(fname_str);
    load( fullfile(folder, subfolder, fname_str) );
    
    % load simulation data
    folder_sim_data = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\data';
    fname_str = sprintf('stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d',...
        num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    fname = fullfile(folder_sim_data, fname_str);
    %disp(fname_str);
    load(fname);
    %}
    %--------------------------------------------------------------------------
    %% Plot phase diagram
    Con_wave = squeeze(Con_all(network, squeeze(trav_wave_conds_met(network, :)), :));
    K_wave = squeeze(K_all(network, squeeze(trav_wave_conds_met(network, :)), :, :));
    
    idx_true = find(trav_wave_all);
    idx_false = find(~trav_wave_all);
    Con_wave_true = Con_wave(idx_true, :);
    Con_wave_false = Con_wave(idx_false, :);
    K_wave_true = K_wave(idx_true, :, :);
    K_wave_false = K_wave(idx_false, :, :);
    
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
            scatter(K_wave_true(:,idx_i,idx_j), Con_wave_true(:,idx_j), 'o',...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.1 0.1 0.1], 'MarkerFaceAlpha', 0.5);
            scatter(K_wave_false(:,idx_i,idx_j), Con_wave_false(:,idx_j), 'o',...
                'MarkerEdgeColor', 'w', 'MarkerFaceColor', [0.9 0.6 0.9], 'MarkerFaceAlpha', 0.5);
            %scatter([500 510 520 530], [100 90 80 110], 'x', 'MarkerEdgeColor', [0.9 0.6 0.9]);
            
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
