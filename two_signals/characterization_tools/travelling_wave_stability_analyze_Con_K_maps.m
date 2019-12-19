% Analyze parameter sets for which TWs can propagate in more detail
% Plot K-Con plots for each interaction, compare 
clear all
close all
set(0, 'defaulttextinterpreter', 'tex')

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

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% specify wave type and characteristics
wave_types_str = {'straight', 'inward bend', 'outward bend'};
wave_type = 1;
num_waves = 1; % number of waves
bandwidth = 1; % width of band of cells of one type

% Data folder
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
%folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
%folder = 'L:\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
%subfolder = 'run1_no_Con_K_info';
subfolder = 'run2';
%subfolder = 'run3_vary_a0_lambda12';
%subfolder = 'run2b_n_pset_10e6';

% save figure folder
%save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_stability';
%save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\temp';
save_folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\all_networks_analytical_phases\phase_diagram';

%% Load analyzed data on TW propagation (analytical)
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
%% Load and process detailed data per network
%for idx_loop=1 %1:numel(x_found)
    idx_loop = 1;
    wave_idx = x_found(idx_loop);
    network = y_found(idx_loop);
    fprintf('wave_idx %d, network %d \n', wave_idx, network);
    states_perm = P(wave_idx, :);
    %}
    %--------------------------------------------------------------------------
    % Get K, Con parameters for waves (analytical)
    % load data
    %folder = 'L:\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general';
    %subfolder = 'run3_vary_a0_lambda12';
    fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
            num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    load( fullfile(folder, subfolder, fname_str) );
    
    Con_all_network = squeeze(Con_all(network, :, :));
    K_all_network = squeeze(K_all(network, :, :, :));
    trav_wave_conds_met_theory = squeeze(trav_wave_conds_met(network, :, :))';
    
    Con_wave_theory = squeeze(Con_all(network, trav_wave_conds_met_theory, :));
    K_wave_theory = squeeze(K_all(network, trav_wave_conds_met_theory, :, :));
    num_wave_sims = size(Con_wave_theory, 1);
    
    clear('Con_all', 'K_all', 'trav_wave_conds_met')
    %M_int = M_int_found{idx_loop}; 
    %%
    % Get K, Con parameters for waves (simulations)
    % load data
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2_stability_sim';
    %if remote
    %    load_folder = strrep(load_folder, 'N:\', 'W:\staff-bulk\');
    %end
    fname_str = sprintf('stability_sim_from_pred_trav_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d',...
           num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));       
    fname = fullfile(load_folder, fname_str);
    load(fname, 'trav_wave_all', 'trav_wave_all_strict', 'unmodified_all');
    
    % get list of simulations that generated trav. waves
    trav_wave_conds_met_sim = (trav_wave_all & unmodified_all);
    fprintf('Network %d, #trav. wave conditions = %d \n', network, sum(trav_wave_conds_met_sim) )
    
    Con_wave_sim = squeeze( Con_all_network(trav_wave_conds_met_sim, :) );
    K_wave_sim = squeeze(K_all_network(trav_wave_conds_met_sim, :, :));    
    
    % get K, Con parameters found by both theory and simulations
    idx_sel = squeeze(trav_wave_conds_met_theory & trav_wave_conds_met_sim);
    Con_wave_both = Con_all_network( idx_sel, :);
    K_wave_both = K_all_network(idx_sel, :, :);  
    
    % get K, Con parameters false positives (theory: TW, simulation: not
    % TW)
    idx_sel = squeeze(trav_wave_conds_met_theory & ~trav_wave_conds_met_sim);
    Con_wave_FP = Con_all_network( idx_sel, :);
    K_wave_FP = K_all_network(idx_sel, :, :);  
    
    % get K, Con parameters false negatives (theory: Not TW, simulation:
    % TW)
    idx_sel = squeeze(~trav_wave_conds_met_theory & trav_wave_conds_met_sim);
    Con_wave_FN = Con_all_network( idx_sel, :);
    K_wave_FN = K_all_network(idx_sel, :, :);  
    
%end
%% Plot specific conditions
%------------Preliminary calculations -------------------------------------
% Parameters
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
%mcsteps = 0;
%[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% -----------Calculate interaction strengths-------------------------------
% calculate fN
[fN, fnn, fnnn] = calc_fN(a0, rcell, gz, dist, lambda);
%{
Rcell = a0*rcell;
idx = gz + round(gz/2); % pick cell not at corner of grid
dist_vec = a0*dist(idx,:);
r = dist_vec(dist_vec>0);
fN = zeros(2,1);
fN(1) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
fN(2) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));

% calculate fnn
fnn = zeros(1, 2);
rnn = a0;
fnn(1) = sinh(Rcell)*exp((Rcell-rnn)./lambda(1)).*(lambda(1)./rnn) ; % calculate signaling strength
fnn(2) = sinh(Rcell)*exp((Rcell-rnn)./lambda(2)).*(lambda(2)./rnn) ; % calculate signaling strength

fnnn = zeros(2,2);
rnnn = [sqrt(3) 2].*a0;
fnnn(1,:) = sinh(Rcell)*exp((Rcell-rnnn)./lambda(1)).*(lambda(1)./rnnn);
fnnn(2,:) = sinh(Rcell)*exp((Rcell-rnnn)./lambda(2)).*(lambda(2)./rnnn);
%}
%------------------------------------ Calculate Y_all----------------------
default_states = [0 0; 0 1; 1 0; 1 1];

% calculate Y(alpha)
% neighbour data
switch wave_type
    case 1
        n_nei = [2	0	0	4; % EF
            2	2	0	2; % F
            2	2	2	0; % M 
            0	2	2	2; % B
            0	0	2	4; % EB
            0	0	0	6]; % E
    case 2
        n_nei = [3	0	0	3;
            2	3	0	1;
            1	2	3	0;
            0	1	2	3;
            0	0	1	5;
            0	0	0	6];
    case 3
        n_nei = [1	0	0	5;
            2	1	0	3;
            3	2	1	0;
            0	3	2	1;
            0	0	3	3;
            0	0	0	6];
end
states = default_states(states_perm, :); % F, M, B, E
types = [states(4,:); states(1,:); states(2,:); states(3,:); states(4,:); states(4,:)];
Y_self = (Con-1).*types + 1;

% calculate Y_nei
Y_nei = zeros(size(n_nei, 1), 2);
for i=1:6
    Y_nei(i,1) = fnn(1)*n_nei(i,:)*((Con(1)-1)*states(:,1)+1);
    Y_nei(i,2) = fnn(2)*n_nei(i,:)*((Con(2)-1)*states(:,2)+1);
end
%
% estimate p
tmp = num_waves*bandwidth;
tmp2 = [tmp tmp tmp gz-3*tmp];
p = (tmp2*states)/gz;

% calculate Y_mf (mean-field contribution)
z = 6; % coordination number
Y_mf = zeros(1, 2);
Y_mf(1) = (fN(1) - z*fnn(1))*( Con(1)*p(1) + (1-p(1)) );
Y_mf(2) = (fN(2) - z*fnn(2))*( Con(2)*p(2) + (1-p(2)) );
%disp(Y_mf)

Y_all = Y_self + Y_nei + Y_mf;

%---------------------------Calculate p values of TWs ---------------------
% Load initial conditions
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_snapshots';
%fname_str = 'trav_wave_single_horizontal_inward_bend';
%fname_str = 'trav_wave_single_horizontal_outward_bend'; 
fname_str = 'trav_wave_single_vertical'; 

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

pw = mean(cells_in, 1);
%--------------------------------------------------------------------------
%% case study: network 15 -> Plot derived bounds only

gz = 50;
mcsteps = 0;
[~, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

a0 = 0.15;
[fN, fnn, ~] = calc_fN(a0, rcell, gz, dist, lambda);

pw = [1/8 1/8];

Con1_all = linspace(1, 1000, 100);
Con2_all = linspace(1, 1000, 100);
YMF1_all = (fN(1) - 6*fnn(1))*(Con1_all*pw(1) + (1-pw(1)));
YMF2_all = (fN(2) - 6*fnn(2))*(Con2_all*pw(2) + (1-pw(2)));

label = 'bounds_only';

% bounds
K11_all_upb = 1 + 2*Con1_all*fnn(1) + 4*fnn(1) + YMF1_all;
K11_all_lwb = 1 + 6*fnn(1) + YMF1_all;
K12_all_upb = Con2_all + (4*Con2_all + 2)*fnn(2) + YMF2_all;
K12_all_lwb = 1 + (2*Con2_all + 4)*fnn(2) + YMF2_all;
K21_all_upb = Con1_all + 4*Con1_all*fnn(1) + 2*fnn(1) + 4*fnn(1) + YMF1_all;
K21_all_lwb = 1 + 2*Con1_all*fnn(1) + 4*fnn(1) + YMF1_all;

% intersection points 
C1_self = (L - 1 - 4.*fnn - (fN' - 6*fnn).*(1-pw))./(2.*fnn+(fN'-6.*fnn).*pw);
C2_self = (L - 1 - 6.*fnn - (fN' - 6*fnn).*(1-pw))./((fN'-6.*fnn).*pw);

C1_mutual = (L - 0 - 2.*fnn - (fN' - 6*fnn).*(1-pw))./(1 + 4.*fnn+(fN'-6.*fnn).*pw);
C2_mutual = (L - 1 - 6.*fnn - (fN' - 6*fnn).*(1-pw))./(2*fnn + (fN'-6.*fnn).*pw);

figure;
subplot(2,2,1);
hold on
plot(K11_all_upb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
plot(K11_all_lwb, Con2_all, '--', 'Color', 'r', 'LineWidth', 2 )
plot([1 L], [C1_self(1) C1_self(1)], '--', 'Color', 'b');
plot([1 L], [C2_self(1) C2_self(1)], '--', 'Color', 'b');
xlim([0 1000]);
ylim([0 1000]);
title('1\leftarrow 1');

subplot(2,2,2);
hold on
plot(K12_all_upb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
plot(K12_all_lwb, Con2_all, '--', 'Color', 'r', 'LineWidth', 2 )
plot([1 L], [C1_mutual(2) C1_mutual(2)], '--', 'Color', 'b');
plot([1 L], [C2_mutual(2) C2_mutual(2)], '--', 'Color', 'b');
xlim([0 1000]);
ylim([0 1000]);
title('1 \leftarrow 2');

subplot(2,2,3);
hold on
plot(K21_all_upb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
plot(K21_all_lwb, Con2_all, '--', 'Color', 'r', 'LineWidth', 2 )
plot([1 L], [C1_mutual(1) C1_mutual(1)], '--', 'Color', 'b');
plot([1 L], [C2_mutual(1) C2_mutual(1)], '--', 'Color', 'b');
xlim([0 1000]);
ylim([0 1000]);
title('2 \leftarrow 1');

%% Find intersections with region boundaries


%% case study: network 15 -> Plot data together with derived bounds
Con1_all = linspace(1, 1000, 100);
Con2_all = linspace(1, 1000, 100);
YMF1_all = (fN(1) - 6*fnn(1))*((Con1_all - 1)*pw(1) + (1-pw(1)));
YMF2_all = (fN(2) - 6*fnn(2))*((Con2_all - 1)*pw(2) + (1-pw(2)));
%label = 'false_and_true_positives'; 
label = 'TP_FP_FN'; 
%label = 'bounds_only';
K_wave_selected = K_wave_both;
Con_wave_selected = Con_wave_both;

% dialogue 15:
% Bounds on K21
% 1 + 2*Con1_all*fnn(1) + 4*fnn(1) + (fN(1) - 6*fnn(1))*((Con1_all - 1)*pw(1) + (1-pw(1))) < K21
% Con1_all + 4*Con1_all*fnn(1) + 2*fnn(1) + 4*fnn(1) + (fN(1) - 6*fnn(1))*(Con1_all - 1)*pw(1) + (1-pw(1)) > K21
% (for 15*, flip 1 and 2)
K21_all_upb = 1 + 2*Con1_all*fnn(1) + 4*fnn(1) + YMF1_all;
K21_all_lwb = Con1_all + 4*Con1_all*fnn(1) + 2*fnn(1) + 4*fnn(1) + YMF1_all;

idx_i = 2; idx_j = 1;
%h = plot_phase_diagram_local(a0, rcell, lambda, dist, idx_i, idx_j);
h = plot_Con_K_map_minimal(idx_i, idx_j);
hold on
box on
%
p1 = scatter(K_wave_selected(:,idx_i,idx_j), Con_wave_selected(:,idx_j), 'k', 'o', 'filled');
p1.MarkerFaceAlpha = 0.5;
p2 = scatter(K_wave_FP(:,idx_i,idx_j), Con_wave_FP(:,idx_j), 'o', 'filled', 'MarkerFaceColor', 'm');
p2.MarkerFaceAlpha = 0.5;
p3 = scatter(K_wave_FN(:,idx_i,idx_j), Con_wave_FN(:,idx_j), 'o', 'filled', 'MarkerFaceColor', 'b');
p3.MarkerFaceAlpha = 0.5;
%}
plot(K21_all_upb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
plot(K21_all_lwb, Con2_all, '-', 'Color', 'r', 'LineWidth', 2 )
%plot([0 1000], [C1(1) C1(1)], '--', 'Color', 'k', 'LineWidth', 2 );
legend([p1 p2 p3], 'True positives', 'False positives', 'False negatives', 'Location', 'se');

xlim([0 1000]);
ylim([0 1000]);

% Save figure
qsave = 0;
fname_str = sprintf('K_Con_plot_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d_interaction_%d_%d_v1_%s',...
	num_waves, wave_type, network, states_perm(1), states_perm(2),...
	states_perm(3), states_perm(4), idx_i, idx_j, label);
fname = fullfile(save_folder, 'with_analytical_bounds', fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);
%%
% Bounds on K12
% 1 + Con2_all*fnn(2) + 4*fnn(2) + (fN(2) - 6*fnn(2))*(Con(2) - 1)*pw(2) + (1-pw(2)) < K12
% Con2_all + 4*Con2_all*fnn(1) + 2*fnn(1) + (fN(2) - 6*fnn(2))*(Con(2) - 1)*pw(2) + (1-pw(2)) > K12

K12_all_lwb = 1 + (2*Con2_all + 4)*fnn(2) + YMF2_all;
K12_all_upb = Con2_all + (4*Con2_all + 2)*fnn(2) + YMF2_all;

idx_i = 1; idx_j = 2;
%h = plot_phase_diagram_local(a0, rcell, lambda, dist, idx_i, idx_j);
h = plot_Con_K_map_minimal(idx_i, idx_j);
hold on
box on
p1 = scatter(K_wave_selected(:,idx_i,idx_j), Con_wave_selected(:,idx_j), 'k', 'o', 'filled');
p1.MarkerFaceAlpha = 0.5;
p2 = scatter(K_wave_FP(:,idx_i,idx_j), Con_wave_FP(:,idx_j), 'o', 'filled', 'MarkerFaceColor', 'm');
p2.MarkerFaceAlpha = 0.5;
p3 = scatter(K_wave_FN(:,idx_i,idx_j), Con_wave_FN(:,idx_j), 'o', 'filled', 'MarkerFaceColor', 'b');
p3.MarkerFaceAlpha = 0.5;

plot(K12_all_upb, Con1_all, '-', 'Color', 'r', 'LineWidth', 2 )
plot(K12_all_lwb, Con1_all, '-', 'Color', 'r', 'LineWidth', 2 )
xlim([0 1000]);
ylim([0 1000]);
legend([p1 p2 p3], 'True positives', 'False positives', 'False negatives', 'Location', 'se');

% Save figure
qsave = 0;
fname_str = sprintf('K_Con_plot_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d_interaction_%d_%d_v1_%s',...
	num_waves, wave_type, network, states_perm(1), states_perm(2),...
	states_perm(3), states_perm(4), idx_i, idx_j, label);
fname = fullfile(save_folder, 'with_analytical_bounds', fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%%
% Bounds on K11
% 1 + 2*Con1_all*fnn(1) + 4*fnn(1) + (fN(1) - 6*fnn(1))*(Con(1) - 1)*pw(1) + (1-pw(1)) > K11
% 1 + 6*fnn(1) + (fN(1) - 6*fnn(1))*(Con(1) - 1)*pw(1) + (1-pw(1)) < K11

K11_all_upb = 1 + 2*Con1_all*fnn(1) + 4*fnn(1) + YMF1_all;
K11_all_lwb = 1 + 6*fnn(1) + YMF1_all;
%K11_all_lwb2 = Con1_all + 4*Con1_all*fnn(1) + 2*fnn(1) + YMF1_all;

idx_i = 1; idx_j = 1;
%h = plot_phase_diagram_local(a0, rcell, lambda, dist, idx_i, idx_j);
h = plot_Con_K_map_minimal(idx_i, idx_j);
hold on
box on


% plot data points
p1 = scatter(K_wave_selected(:,idx_i,idx_j), Con_wave_selected(:,idx_j), 'k', 'o', 'filled');
p1.MarkerFaceAlpha = 0.5;
p2 = scatter(K_wave_FP(:,idx_i,idx_j), Con_wave_FP(:,idx_j), 'o', 'filled', 'MarkerFaceColor', 'm');
p2.MarkerFaceAlpha = 0.5;
p3 = scatter(K_wave_FN(:,idx_i,idx_j), Con_wave_FN(:,idx_j), 'o', 'filled', 'MarkerFaceColor', 'b');
p3.MarkerFaceAlpha = 0.5;

% plot bounds
p11 = plot(K11_all_upb, Con1_all, '-', 'Color', 'r', 'LineWidth', 2 );
p12 = plot(K11_all_lwb, Con1_all, '-', 'Color', 'r', 'LineWidth', 2 );
legend([p1 p2 p3], {'True positives', 'False positives', 'False negatives'}, 'Location', 'se');

%plot(K11_all_lwb2, Con1_all, 'c--', 'LineWidth', 2 )
%legend([p1 p], 'False positives', 'True positives', 'Location', 'se');
xlim([0 1000]);
ylim([0 1000]);

% Save figure
qsave = 0;
fname_str = sprintf('K_Con_plot_wave_num_%d_type_%d_network_%d_states_F%d_M%d_B%d_E%d_interaction_%d_%d_v1_%s',...
	num_waves, wave_type, network, states_perm(1), states_perm(2),...
	states_perm(3), states_perm(4), idx_i, idx_j, label);
fname = fullfile(save_folder, 'with_analytical_bounds', fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot phase space regions bounds (analytical)

% calculate fN, fnn
%{
this_a0 = 1.5;
gz_all = 5:5:50;
fN_all = zeros(numel(gz_all), 2);
for i=1:numel(gz_all)
    this_gz = gz_all(i);
    [~, this_dist] = initial_cells_random_markov_periodic(this_gz, mcsteps, rcell);
    [fN, fnn, ~] = calc_fN(this_a0, rcell, this_gz, this_dist, lambda);
    fN_all(i, :) = fN;
end

h=figure;
plot(gz_all, fN_all);
%}

%% Check conditions for specific parameter set

% to delete
%{
conditions_met = zeros(6,1);
% (1) E_F->F
% (2) F->M
% (3) M->B
% (4) B->E
% (5) E_B->E
% (6) E->E 
targets = [1 2 3 4 4 4]; % recall F, M, B, E
output_state_list = zeros(6, 2);
for i=1:6
    %i = 1;
    Y = Y_all(i, :); % 2x1
    
    %output_state = prod(((Y - K).*M_int>0) + (1-abs(M_int)), 2)';
    target_state = states(targets(i),:);
    conditions_met(i) = all(output_state==target_state);
    output_state_list(i,:) = output_state;
end
%disp(output_state_list);
%disp(conditions_met);

trav_wave_conds_met = all(conditions_met);
    
%}
 %% calculate fN
    Rcell = a0*rcell;
    idx = gz + round(gz/2); % pick cell not at corner of grid
    dist_vec = a0*dist(idx,:);
    r = dist_vec(dist_vec>0);
    fN = zeros(2,1);
    fN(1) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
    fN(2) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));

    % calculate fnn
    fnn = zeros(1, 2);
    rnn = a0;
    fnn(1) = sinh(Rcell)*exp((Rcell-rnn)./lambda(1)).*(lambda(1)./rnn) ; % calculate signaling strength
    fnn(2) = sinh(Rcell)*exp((Rcell-rnn)./lambda(2)).*(lambda(2)./rnn) ; % calculate signaling strength

    fnnn = zeros(2,2);
    rnnn = [sqrt(3) 2].*a0;
    fnnn(1,:) = sinh(Rcell)*exp((Rcell-rnnn)./lambda(1)).*(lambda(1)./rnnn);
    fnnn(2,:) = sinh(Rcell)*exp((Rcell-rnnn)./lambda(2)).*(lambda(2)./rnnn);
%% Functions

function [fN, fnn, fnnn] = calc_fN(a0, rcell, gz, dist, lambda)
    % calculate fN
    Rcell = a0*rcell;
    idx = gz + round(gz/2); % pick cell not at corner of grid
    dist_vec = a0*dist(idx,:);
    r = dist_vec(dist_vec>0);
    fN = zeros(2,1);
    fN(1) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
    fN(2) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));

    % calculate fnn
    fnn = zeros(1, 2);
    rnn = a0;
    fnn(1) = sinh(Rcell)*exp((Rcell-rnn)./lambda(1)).*(lambda(1)./rnn) ; % calculate signaling strength
    fnn(2) = sinh(Rcell)*exp((Rcell-rnn)./lambda(2)).*(lambda(2)./rnn) ; % calculate signaling strength

    fnnn = zeros(2,2);
    rnnn = [sqrt(3); 2].*a0;
    fnnn(1,:) = sinh(Rcell)*exp((Rcell-rnnn)./lambda(1)).*(lambda(1)./rnnn);
    fnnn(2,:) = sinh(Rcell)*exp((Rcell-rnnn)./lambda(2)).*(lambda(2)./rnnn);
end
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

function h = plot_Con_K_map_minimal(idx_gene, idx_mol)
    % Plot phase diagram without colors
    % idx_gene: index of gene under control
    % idx_mol: index of sensed molecule / molecule affecting the gene (1 <= i <= L)
        
    % Parameter range map
    Con_max = 1000;
    K_max = 1000;

    % plot figure
    h = figure;
    hold on
    %set(him, 'AlphaData', out > 0); % invisible if not from any region
    
    % adjust the graph
    set(gca,'ydir', 'normal', 'FontSize', 24)
    xlabel(sprintf('K^{(%d%d)}', idx_gene, idx_mol), 'FontSize', 24)
    ylabel(sprintf('C_{ON}^{(%d)}', idx_mol), 'FontSize', 24)
    ylim([1 Con_max])
    xlim([1 K_max])
    %title(sprintf('$$f_N = %.3f, a_0 = %.2f$$', fN, a0), 'FontSize', 30)
    title(sprintf('Interaction %d \\leftarrow %d', idx_gene, idx_mol))
    %title(sprintf('Molecule %d', idx_mol));    
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
    ylabel('$$C_{ON}$$', 'FontSize', 24)
    ylim([1 Con_max])
    xlim([1 K_max])
    %title(sprintf('$$f_N = %.3f, a_0 = %.2f$$', fN, a0), 'FontSize', 30)
    title(sprintf('Interaction $$%d \\leftarrow %d$$', idx_gene, idx_mol))
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