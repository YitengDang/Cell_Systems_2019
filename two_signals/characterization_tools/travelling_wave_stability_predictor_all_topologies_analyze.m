% Predicts whether a travelling wave can propagate according to nearest
% neighbour interactions
% Apply computational search over all topologies to try to find parameters
% that allow for propagation of any kind of travelling wave
clear all
close all
set(0,'defaulttextinterpreter', 'latex')

%% Parameters
% Number of parameter sets to do
n_pset = 10^6;

% Manual input
gz = 15;
N = gz^2;
%a0 = 1.5;
rcell = 0.2;
%Rcell = rcell*a0;
%lambda = [1 1.2];
%M_int = [0 1; -1 1]; % network 15 reversed
%{
M_int = [0 1; -1 1]; % network 15 reversed
M_int = [1 -1; 1 0]; % network 15
M_int = [1 1; -1 0]; % network 19
M_int = [1 -1; 1 1]; % network 33
M_int = [-1 -1; 1 1]; % network 34
M_int = [-1 1; -1 1]; % network 36
%}
%Con = [18 16];
%K = [0 9; 11 4];
%K = [4 11; 9 0];

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Obtain from simulation
%{
folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected';
subfolder = 'patterns\Network 33';
fname_str = 'tripple_wave_diagonally';
fname = fullfile(folder, subfolder, fname_str);
load(fname, 'save_consts_struct');

Con = save_consts_struct.Con;
K = save_consts_struct.K;
N = save_consts_struct.N;
gz = sqrt(N);
a0 = save_consts_struct.a0;
rcell = save_consts_struct.rcell;
Rcell = rcell*a0;
lambda = save_consts_struct.lambda;
M_int = save_consts_struct.M_int;
%}

% specify wave type and characteristics
wave_types_str = {'straight', 'inward bend', 'outward bend'};
wave_type = 1;
num_waves = 1; % number of waves
bandwidth = 1; % width of band of cells of one type

% order of states: F, M, B, E
states_perm = [1 2 3 4];
%{
states_perm = [2 4 3 1]; % network 15 reversed
states_perm = [3 4 2 1]; % network 15 
states_perm = [4 3 1 2]; % network 19
states_perm = [3 4 2 1]; % network 33
states_perm = [4 2 1 3]; % network 33/34
states_perm = [2 4 3 1]; % network 36
default_states = [0 0; 0 1; 1 0; 1 1];
%}

% Data folder
folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
%folder = 'L:\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
%subfolder = 'run1_no_Con_K_info';
%subfolder = 'run2';
subfolder = 'run3_vary_a0_lambda12';
    
%% 
%{
types_waves = 24;
num_networks = 44;
wave_possible = zeros(types_waves, num_networks); % stores whether a wave of a certain type is possible in a given network
P = perms(1:4);

for i=1:types_waves
    disp(i);
    states_perm = P(i, :);
    
    % Load data
    
    fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
            num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    load( fullfile(folder, subfolder, fname_str) );
    
    % Find networks that can support the wave
    % (1) For old set (old network numbering 1-81)
    %{
    networks_idx_rev = zeros(81, 1);
    networks_idx_rev(networks_idx) = 1:44;
    idx_network_1 = find( sum(trav_wave_conds_met, 2)>0 ); % network indices, 1:81 indexing
    idx_network_2 = networks_idx_rev( idx_network_1 );
    %}
    
    % (2) For new set 
    idx_network_2 = find( sum(trav_wave_conds_met, 2)>0 ); % networks, new indexing
    disp('Found networks:');
    disp(idx_network_2);
    disp('----------------');
    
    % To do
    if ~isempty(idx_network_2)
        wave_possible(i, idx_network_2) = 1;
    end
end
%}
%% Save analyzed data
%{
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_%s',...
    num_waves, wave_type, subfolder);
save( fullfile(folder, fname_str) );
%}
%% Load analyzed data
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_%s',...
    num_waves, wave_type, subfolder);
load( fullfile(folder, fname_str) );

%% Plot wave_possible as heatmap
h = figure;
imagesc(1:44, 1:24, wave_possible);
set(gca, 'YDir', 'normal', 'FontSize', 20);
xlabel('Network');
ylabel('Wave type');

qsave = 1;
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_imagesc_%s',...
    num_waves, wave_type, subfolder);
fname = fullfile(folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% display found waves
[x_found, y_found] = find(wave_possible);
t=table(P(x_found, :), y_found, 'VariableNames', {'Wave_type', 'Network'});
t2=table(x_found, y_found, 'VariableNames', {'wave_idx', 'Network'});

disp('Cell states F, M, B, E');
disp(t);
disp(t2);

%% Find properties of specific travelling wave
Qvals = zeros(numel(x_found), 1); % Q-values of circuits
for idx_loop=1:numel(x_found)
    %wave_idx = 2;
    %network = 19;
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    
    states_perm = P(wave_idx, :);
    %
    % load data
    folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
    %subfolder = 'run3_vary_a0_lambda12';
    fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
            num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    load( fullfile(folder, subfolder, fname_str) );
    %}
    %% Number / fraction of parameter sets giving waves
    %
    %n_pset_wave = sum(trav_wave_conds_met(networks_idx(network), :));
    n_pset_wave = sum(trav_wave_conds_met(network, :));
    frac_pset_wave = n_pset_wave/n_pset;
    fprintf('Fraction of parameter sets yielding waves: %.10f \n', frac_pset_wave);
    Qvals(idx_loop) = n_pset_wave; % store absolute numbers rather than fractions
    %}
    %% parameter values of simulations giving waves
    Con_wave = squeeze(Con_all(network, squeeze(trav_wave_conds_met(network, :)), :));
    K_wave = squeeze(K_all(network, squeeze(trav_wave_conds_met(network, :)), :, :));

    % Find structure in pset
    h = figure;
    hold on
    histogram(Con_wave(:,1), 1:50:1000, 'normalization', 'pdf');
    histogram(Con_wave(:,2), 1:50:1000, 'normalization', 'pdf');
    legend({'C_{ON}^{(1)}', 'C_{ON}^{(2)}'}, 'Location', 'best');
    xlabel('Value');
    ylabel('Probability');
    title('$C_{ON}$ values of travelling waves');
    set(gca, 'FontSize', 20);

    h2 = figure;
    hold on
    histogram(K_wave(:,1, 1), 1:50:1000, 'normalization', 'pdf');
    histogram(K_wave(:,1, 2), 1:50:1000, 'normalization', 'pdf');
    histogram(K_wave(:,2, 1), 1:50:1000, 'normalization', 'pdf');
    histogram(K_wave(:,2, 2), 1:50:1000, 'normalization', 'pdf');
    legend({'K^{(1,1)}', 'K^{(1,2)}', 'K^{(2,1)}', 'K^{(2,2)}'});
    xlabel('Value');
    ylabel('Probability');
    title('$K$ values of travelling waves');
    set(gca, 'FontSize', 20);
    
    %
    % save figures
    qsave = 1;
    fname_str = sprintf('Trav_wave_conditions_wave_num_%d_type_%d_networks_%d_states_F%d_M%d_B%d_E%d',...
        num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    fname = fullfile(folder, strcat(fname_str, '_Con_hist'));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
    fname = fullfile(folder, strcat(fname_str, '_K_hist'));
    save_figure(h2, 10, 8, fname, '.pdf', qsave);
    
    %% Calculate correlations
    if idx_loop<3 % networks with 3 links
        A = [Con_wave K_wave(:,1:3)]; % variables: Con(1) Con(2) K(1,1) K(2,1) K(1,2) K(2,2)
        R = corrcoef(A);

        h3=figure;
        imagesc(R);
        plot_labels = {'C_{ON}^{(1)}', 'C_{ON}^{(2)}',...
            'K^{(1,1)}', 'K^{(2,1)}', 'K^{(1,2)}'};
        set(gca, 'XTick', 1:5, 'XTickLabels', plot_labels);
        set(gca, 'YTick', 1:5, 'YTickLabels', plot_labels);
    else % networks with 4 links 
        A = [Con_wave K_wave(:,:)]; % variables: Con(1) Con(2) K(1,1) K(2,1) K(1,2) K(2,2)
        R = corrcoef(A);

        h3=figure;
        imagesc(R);
        plot_labels = {'C_{ON}^{(1)}', 'C_{ON}^{(2)}',...
            'K^{(1,1)}', 'K^{(2,1)}', 'K^{(1,2)}', 'K^{(2,2)}'};
        set(gca, 'XTick', 1:6, 'XTickLabels', plot_labels);
        set(gca, 'YTick', 1:6, 'YTickLabels', plot_labels);
    end
    caxis([-1 1]);
    set(gca, 'FontSize', 20);
    
    % make a new colormap
    c = colorbar;
    ylabel(c, '\rho');
    greenColorMap = [zeros(1, 130), linspace(0, 1, 126).^(1/2)];
    redColorMap = [linspace(1, 0, 126).^(1/2), zeros(1, 130)];
    colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
    colormap(colorMap);
    
    % Save correlation heatmap
    qsave = 1;
    fname_str = sprintf('Trav_wave_conditions_wave_num_%d_type_%d_networks_%d_states_F%d_M%d_B%d_E%d',...
        num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    fname = fullfile(folder, strcat(fname_str, '_Con_K_corr'));
    save_figure(h3, 10, 8, fname, '.pdf', qsave);
end

%% Plot Q values as bar graph
h = figure;
labels = {'15, [3 4 2 1]','19, [4 3 1 2]','33, [4 2 1 3]','33, [3 4 2 1]',...
    '34, [4 2 1 3]','36, 2 4 3 1]'};
c = categorical(labels);
bar(c, Qvals/n_pset);
xlabel('Network, wave form')
ylabel('Q-value');
set(gca, 'FontSize', 20);

% save figures
qsave = 0;
folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
fname_str = sprintf('Trav_wave_conditions_wave_num_%d_type_%d',...
    num_waves, wave_type);
fname = fullfile(folder, strcat(fname_str, '_Q_vals_by_wave'));
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Analyze waves of a specific wave_idx and network
% Todo: *quantify by calculating fraction in each phase region
% *Look at combined phases for all interactions together
% *Refer back to state diagrams: do these predict the possible phases?

for idx_loop=2:6
    %idx_loop = 1;
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    states_perm = P(wave_idx, :);

    % load data
    folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
    %subfolder = 'run3_vary_a0_lambda12';
    fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
            num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    load( fullfile(folder, subfolder, fname_str) );
    %% Plot phase diagram
    Con_wave = squeeze(Con_all(network, squeeze(trav_wave_conds_met(network, :)), :));
    K_wave = squeeze(K_all(network, squeeze(trav_wave_conds_met(network, :)), :, :));
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
            %hold on
            scatter(Con_wave(:,idx_j), K_wave(:,idx_i,idx_j), 'k', 'x');

            % save figures
            qsave = 1;
            folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\phase_diagram';
            fname_str = sprintf('Trav_wave_conditions_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d_interaction_%d_%d_scatter',...
                num_waves, wave_type, network, states_perm(1), states_perm(2),...
                states_perm(3), states_perm(4), idx_i, idx_j);
            fname = fullfile(folder, fname_str);
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