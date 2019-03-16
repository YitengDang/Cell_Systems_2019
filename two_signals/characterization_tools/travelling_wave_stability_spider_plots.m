% Tests whether a travelling wave can propagate by performing explicit
% simulations
clear all
close all
set(0, 'defaulttextinterpreter', 'tex');

%% Parameters
% remote destination (Webdrive)?
remote = 0;

% parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];

%Con = [18 16];
%Coff = [1 1];
%M_int = [0 1; -1 1];
%K = [0 9; 11 6];

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
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_snapshots';
if remote
    load_folder = strrep(load_folder, 'N:\', 'W:\staff-bulk\');
end
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
%% Load overview data
% Load list of found networks and wave forms
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
if remote
    folder = strrep(folder, 'N:\', 'W:\staff-bulk\');
end
subfolder = 'run2';
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_%s',...
    num_waves, wave_type, subfolder);
load( fullfile(folder, fname_str) );

% get waveforms and networks of possible waves
[x_found, y_found] = find(wave_possible);
num_psets = numel(x_found);
t=table(P(x_found, :), y_found, 'VariableNames', {'Wave_type', 'Network'});
t2=table(x_found, y_found, 'VariableNames', {'wave_idx', 'Network'});

disp('Cell states F, M, B, E');
disp(t);
disp(t2);

%% Plot parameters of simulations giving waves as spider plot
for idx_loop=1:num_psets
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    states_perm = P(wave_idx, :);
    
    % load predictor data to get parameter sets
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2';
    if remote
        load_folder = strrep(load_folder, 'N:\', 'W:\staff-bulk\');
    end
    fname_str = sprintf('Trav_wave_predictor_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d_processed',...
           num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4)); 
    fname = fullfile(load_folder, fname_str);
    
    % ----- Option 1: plot simulation data
    %
    load(fname, 'Con_all', 'K_all', 'trav_wave_conds_met');  % load parameter sets and data
    
    % save name, simulations
    label = 'simulations';
    subfolder = 'spider_plots_simulations';
    %}
    % ----- Option 2: plot analytical data
    %{
    load(fname, 'Con_all', 'K_all'); % load parameter sets
    
    % load simulation data
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2_stability_sim';
    if remote
        load_folder = strrep(load_folder, 'N:\', 'W:\staff-bulk\');
    end
    fname_str = sprintf('stability_sim_from_pred_trav_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d',...
           num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));       
    fname = fullfile(load_folder, fname_str);
    load(fname, 'trav_wave_all', 'trav_wave_all_strict', 'unmodified_all');
    
    % get list of simulations that generated trav. waves
    trav_wave_conds_met = (trav_wave_all & unmodified_all);
    fprintf('Network %d, #trav. wave conditions = %d \n', network, sum(trav_wave_conds_met) )
    
    % save name, analytical
    label = 'analytical'; 
    subfolder = 'spider_plots_analytical';
    %}
    %%
    % filter data on TW simulations
    Con_wave_sim = Con_all(trav_wave_conds_met, :);
    K_wave_sim = K_all(trav_wave_conds_met, :, :);
    
    % Plot data as spider plot
    if idx_loop<3
        K_idx = 1:3;
    else
        K_idx = 1:4;
    end
    
    %P_data = log10([Con_wave_sim, K_wave_sim(:,K_idx)]); % -> more generally, filter on M_int
    P_data = [Con_wave_sim, K_wave_sim(:,K_idx)];
    
    %P_labels = {'$C_{ON}^{(1)}$', '$C_{ON}^{(2)}$', '$K^{(11)}$',...
    %    '$K^{(12)}$', '$K^{(21)}$', '$K^{(22)}$'};
    P_labels = {'C_{ON}^{(1)}', 'C_{ON}^{(2)}', 'K^{(11)}',...
        'K^{(12)}', 'K^{(21)}', 'K^{(22)}'};
    axes_interval = 3;
    
    %{
    if idx_loop==3
        spider_plot(P_data, P_labels([1:2 K_idx+2]), axes_interval,...
            'Marker', 'o',...
            'LineStyle', '-',...
            'Color', [0 0 1],...
            'LineWidth', 2,...
           'MarkerSize', 2);
    else
    %}
        spider_plot_linear(P_data, P_labels([1:2 K_idx+2]), axes_interval,...
            'Marker', 'o',...
            'LineStyle', 'none',...
            'Color', [1 0 0],...
            'LineWidth', 2,...
            'MarkerSize', 2);
    %end
    %
    %title(sprintf('n=%d, nw %d', size(P_data, 1), network),...
    %    'Fontweight', 'bold',...
    %    'FontSize', 28);
    %}
    set(gcf, 'Units', 'Inches', 'Position', [2 2 10 8]);
    h = gcf;
    hold off 
    %hold on
    
    % Save figure
    qsave = 1;
    fname_root = sprintf('Wave_type_%d_network_%d_states_F%d_M%d_B%d_E%d_TW_n%d_%s',...
        1, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4), size(P_data, 1), label);
    save_folder_fig = ...
        fullfile('H:\My Documents\Multicellular automaton\figures\trav_wave_stability', subfolder);
    if remote
        save_folder_fig = strrep(save_folder_fig, 'H:\', 'W:\staff-homes\d\yitengdang\');
    end
    save_figure(h, 10, 8, fullfile(save_folder_fig, ...
        strcat(fname_root, '_TW_params_spider_linear_lines_filled')), '.pdf', qsave);
    
    % close all
    %}
    
    % Save data: filtered set of parameters that give rise to TW propagation (from simulations)
    %{
    folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2_net_parameters_TW_sim';
    if remote
        folder = strrep(folder, 'N:\', 'W:\staff-bulk\');
    end
    fname_str = strcat(fname_root, '_Con_K_values_waves_sim_strict');
    save(fullfile(folder, fname_str), 'Con_wave_sim', 'K_wave_sim', 'N', 'a0', 'rcell', 'lambda', 'hill', 'noise');
    %}
    
    %pause(1);
    %close all;
end

%%
% Save figure
qsave = 0;
save_folder_fig = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\all_networks_test_analytical_in_sims';
fname_root = sprintf('Wave_type_%d_network_%d_both_wave_types',...
    1, network);
save_figure(h, 10, 8, fullfile(save_folder_fig, ...
    strcat(fname_root, '_parameters_wave_prop_spider_layout2')), '.eps', qsave);
% blue -> 
% red -> 
%}
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

%% 
% Earlier loading data
% v1 data
%{
load_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\data';
%fname_str_data = sprintf('%s_stability_sim_%s_K12_K21_K22_range', fname_str, fname_str_default);  
fname_str_data = sprintf('%s_stability_sim_K12_K21_K22_range', fname_str);  
fname = fullfile(load_folder, fname_str_data);
load(fname, 'K_12_all', 'K_21_all', 'K_22_all', 'trav_wave_all', 'unmodified_all');
%}

% v2 data
%{
fname_str_temp = {'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_15_states_3_4_2_1',...
'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_19_states_4_3_1_2',...
'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_33_states_3_4_2_1',...
'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_33_states_4_2_1_3',...
'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_34_states_4_2_1_3',...
'stability_sim_from_pred_trav_wave_single_vertical_wave_num_%d_type_%d_network_36_states_2_4_3_1'};
fname_str_all = cellfun(@(x) sprintf(x, num_waves, wave_type), fname_str_temp, 'UniformOutput', false);
%}
%load_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\data';
%load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
%num_files = numel(fname_str_all);