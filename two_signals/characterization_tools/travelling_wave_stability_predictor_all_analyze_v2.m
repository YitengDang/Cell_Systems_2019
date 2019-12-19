% Predicts whether a travelling wave can propagate according to nearest
% neighbour interactions
% Apply computational search over all topologies to try to find parameters
% that allow for propagation of any kind of travelling wave
% v2: significantly shorten plots, only include relevant ones
clear all
close all
set(0,'defaulttextinterpreter', 'tex')

%% Parameters
% Number of parameter sets to do
n_pset = 10^5;

% Manual input
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];
Coff = [1 1];
logic = 0;
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

% logic str
if logic
    logic_str = 'AND_logic';
else
    logic_str = 'OR_logic';
end

% specify wave type and characteristics
wave_types_str = {'straight', 'inward bend', 'outward bend'};
wave_type = 1;
num_waves = 1; % number of waves
bandwidth = 1; % width of band of cells of one type

% order of states: F, M, B, E
%states_perm = [1 2 3 4];
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
%folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
%folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
%folder = 'L:\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\TW_propagation_conditions_analytical';
%subfolder = 'run1_no_Con_K_info';
%subfolder = 'run2';
%subfolder = 'run3_vary_a0_lambda12';
%subfolder = 'run2b_n_pset_10e6';
subfolder = 'run_stability';
load_folder = fullfile(folder, subfolder);

% save figure folder
%save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_stability';
%save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\temp';
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\TW_propagation_conditions_analytical';
%% Find networks capable of supporting travelling waves 
%
types_waves = 24;
num_networks = 44;
wave_possible = zeros(types_waves, num_networks); % stores whether a wave of a certain type is possible in a given network
P = perms(1:4);

for i=18 %1:types_waves
    disp(i);
    states_perm = P(i, :);
    
    % Load data
    %fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
    % 	num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    fname_str = sprintf('TW_conds_%s_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
        logic_str, num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    load( fullfile(load_folder, fname_str) );
    
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
%
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_%s',...
    num_waves, wave_type, subfolder);
save( fullfile(folder, fname_str) );
%}
%% Load analyzed data
%
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_%s',...
    num_waves, wave_type, subfolder);
load( fullfile(folder, fname_str) );
%}

%% display found waves
%
[x_found, y_found] = find(wave_possible);
t=table(P(x_found, :), y_found, 'VariableNames', {'Wave_type', 'Network'});
t2=table(x_found, y_found, 'VariableNames', {'wave_idx', 'Network'});

disp('Cell states F, M, B, E');
disp(t);
disp(t2);

% Get a list of all interaction matrices for the found networks
[M_int_found] = get_found_M_int(y_found);
%}
%% Plot wave_possible as heatmap
%{
h = figure;
imagesc(1:44, 1:24, wave_possible);
set(gca, 'YDir', 'normal', 'FontSize', 20);
xlabel('Network');
ylabel('Wave type');

qsave = 1;
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_imagesc_%s',...
    num_waves, wave_type, subfolder);
fname = fullfile(save_folder, 'all_networks_analytical_Con_K_conditions', fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%}

%% Calculate Q-values for each wave form and network
Qvals = zeros(numel(x_found), 1); % Q-values of circuits
for idx_loop=1:numel(x_found)
    %wave_idx = 2;
    %network = 19;
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    
    states_perm = P(wave_idx, :);
    %
    % load data - same data files but now only the parameter sets that give waves
    % (1) original AND-logic
    %fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
    % 	num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    % (2) OR-logic
    fname_str = sprintf('TW_conds_%s_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
        logic_str, num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    load( fullfile(load_folder, fname_str) );
    %} 
    %% Number / fraction of parameter sets giving waves
    %
    %n_pset_wave = sum(trav_wave_conds_met(networks_idx(network), :));
    n_pset_wave = sum(trav_wave_conds_met(network, :));
    frac_pset_wave = n_pset_wave/n_pset;
    fprintf('Fraction of parameter sets yielding waves: %.10f \n', frac_pset_wave);
    Qvals(idx_loop) = n_pset_wave; % store absolute numbers rather than fractions
    %}
end

%% Plot Q values as bar graph
h = figure;
hold on
labels = {};
for i=1:numel(x_found)
    labels{i} = sprintf('%d, [%d %d %d %d]', y_found(i), P(x_found(i), 1),...
        P(x_found(i), 2), P(x_found(i), 3), P(x_found(i), 4) );
end
%labels = {'15, [3 4 2 1]','19, [4 3 1 2]','33, [4 2 1 3]','33, [3 4 2 1]',...
%    '34, [4 2 1 3]','36, [2 4 3 1]'};
c = categorical(labels);

%yyaxis left
bar(c, Qvals/n_pset); % unnormalized
ylabel('Robustness'); 
%{
yyaxis right
num_params = [5 5 6 6 6 6];
bar(c, (Qvals/n_pset).^(1./num_params')); % unnormalized
ylabel('Normalized robustness');
%}
xlabel('Network, wave form')
set(gca, 'FontSize', 32);

% save figures
qsave = 0;
%folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
fname_str = sprintf('TW_conds_%s_wave_num_%d_type_%d_%s',...
    logic_str, num_waves, wave_type, subfolder);
fname = fullfile(save_folder,...
    strcat(fname_str, '_Q_vals_by_wave_v2_both'));
save_figure(h, 10, 8, fname, '.pdf', qsave);
%}
%% Plot Con, K values for waves in projections
for idx_loop=1:numel(x_found)
    %idx_loop = 1;
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    states_perm = P(wave_idx, :);
    
    %idx_wave = 1;
    
    % load data
     % (1) original AND-logic
    %fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
    % 	num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    % (2) OR-logic
    fname_str = sprintf('TW_conds_%s_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
        logic_str, num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    load( fullfile(folder, subfolder, fname_str) );
    %}
    %--------------------------------------------------------------------------
    %% Plot parameters in projections
    Con_wave = squeeze(Con_all(network, squeeze(trav_wave_conds_met(network, :)), :));
    K_wave = squeeze(K_all(network, squeeze(trav_wave_conds_met(network, :)), :, :));
    
    h = figure;
    for idx_i=1:2
        for idx_j=1:2
            %idx_i = 1;
            %idx_j = 1;
            
            % check that the interaction is not absent
            if K_wave(1, idx_i, idx_j)==0
                continue
            end
            
            % plot full phase diagram
            %h = plot_phase_diagram_local(a0, rcell, lambda, dist, idx_i, idx_j);
            
            % plot data points only
            %h=figure;
            idx_plot = 2*(idx_i-1) + (idx_j-1) + 1;
            subplot(2, 2, idx_plot);
            scatter(K_wave(:,idx_i,idx_j), Con_wave(:,idx_j), 'k', 'x');
            set(gca, 'FontSize', 32)
            xlabel('K')
            ylabel('C_{ON}')
            ylim([1 1000])
            xlim([1 1000])
            title(sprintf('Interaction %d \\leftarrow %d', idx_i, idx_j))
            
            
        end
    end
    
    % save figures
    qsave = 1;
    %save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\all_networks_analytical_phases\phase_diagram';
    fname_str = sprintf('Trav_wave_conditions_filtered_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d_scatter',...
        num_waves, wave_type, network, states_perm(1), states_perm(2),...
        states_perm(3), states_perm(4));
    fname = fullfile(save_folder, fname_str);
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end
    
%% Functions
function [M_int_found] = get_found_M_int(y_found)
    % Get a list of all interaction matrices for the found networks
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

    M_int_found = cell(numel(y_found), 1);
    for i=1:numel(y_found)
        M_int_found{i} = M_int_all_reduced{y_found(i)};
        %disp(M_int_found{i})
    end
end

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
    ylabel('$C_{ON}$', 'FontSize', 24)
    ylim([1 Con_max])
    xlim([1 K_max])
    %title(sprintf('$$f_N = %.3f, a_0 = %.2f$$', fN, a0), 'FontSize', 30)
    title(sprintf('Interaction $%d \\leftarrow %d$', idx_gene, idx_mol))
    %title(sprintf('Molecule %d', idx_mol));
end

function A = get_state_diagram_from_phase(phase, M_int)
    % essentially second part of plot_state_diagram_multicell

    % Map from phase to diagram
    % activation/repression | state | input molecule (1/2)
    g_map = cell(2, 6, 2);
    % 0=OFF, 1:ON, 2:UNKNOWN
    % activation 
    g_map{1,1,1} = 2*ones(2);
    g_map{1,1,2} = 2*ones(2);
    g_map{1,2,1} = ones(2);
    g_map{1,2,2} = ones(2);
    g_map{1,3,1} = [2 2; 1 1];
    g_map{1,3,2} = [2 1; 2 1];
    g_map{1,4,1} = [0 0; 2 2];
    g_map{1,4,2} = [0 2; 0 2];
    g_map{1,5,1} = zeros(2);
    g_map{1,5,2} = zeros(2);
    g_map{1,6,1} = [0 0; 1 1];
    g_map{1,6,2} = [0 1; 0 1];
    % repression 
    %(note: this is precisely NOT g_map{1,:,:} in the three-val
    % boolean algebra with NOT 2 = 2)
    g_map{2,1,1} = 2*ones(2);
    g_map{2,1,2} = 2*ones(2);
    g_map{2,2,1} = zeros(2);
    g_map{2,2,2} = zeros(2);
    g_map{2,3,1} = [2 2; 0 0];
    g_map{2,3,2} = [2 0; 2 0];
    g_map{2,4,1} = [1 1; 2 2];
    g_map{2,4,2} = [1 2; 1 2];
    g_map{2,5,1} = ones(2);
    g_map{2,5,2} = ones(2);
    g_map{2,6,1} = [1 1; 0 0];
    g_map{2,6,2} = [1 0; 1 0];

    gij = cell(2);
    X_out = cell(2, 1);
    for i=1:2
        if all(M_int(i,:)==0)
            fprintf('No input for gene %d \n', i);
            % no input => output=initial state
            X1_in = [0 0; 1 1]; 
            X2_in = [0 1; 0 1];
            X_in = (i==1).*X1_in + (i==2).*X2_in; 
            X_out{i} = X_in;
        else
            % normal case
            for j=1:2
                if M_int(i,j)~=0
                    idx = (M_int(i,j)==1) + (M_int(i,j)==-1)*2;
                    gij{i,j} = g_map{idx, phase(i,j), j};
                else
                    gij{i,j} = ones(2); % Fixed ambiguous inputs
                end
            end
            X_out{i} = and3(gij{i,1}, gij{i,2}); % three-valued logic
        end
    end

    A = zeros(4); % graph adjacency matrix
    for i=1:2
        for j=1:2
            state_in = i + 2*(j-1); 
            X_out_this = [X_out{1}(i,j) X_out{2}(i,j)]; % tentative 
            %disp(X_out_this)
            if all(X_out_this~=2) % unambiguous out state
                state_out = X_out_this(1)+1 + 2*X_out_this(2); % (i,j) -> idx
                A(state_in, state_out) = 1;
                %disp(state_out);
            elseif sum(X_out_this==2)==1 % semi-definite
                if (X_out_this(1)==2)
                    X_out_both = [0 X_out_this(2); 1 X_out_this(2)];
                elseif (X_out_this(2)==2)
                    X_out_both = [X_out_this(1) 0; X_out_this(1) 1];
                end
                state_out = X_out_both*[1; 2]+1;
                %[X_out_both(1,1)+1 + 2*X_out_both(1,2);...
                %    X_out_both(2,1)+1 + 2*X_out_both(2,2)];
                A(state_in, state_out) = 1;
                %disp(state_out);
            elseif sum(X_out_this==2)==2 
                A(state_in, :) = 1;
            end
        end
    end
end

function out = and3(x,y)
    out = min(x.*y, 2);
end