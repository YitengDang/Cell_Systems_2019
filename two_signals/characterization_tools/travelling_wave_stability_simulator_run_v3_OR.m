% Tests whether a travelling wave can propagate by performing explicit
% simulations
% v3: run for ALL parameter sets that the predictor tested
% Needed to obtain precision and recall scores
% OR: adapted for OR logic gate
clear all
close all
set(0,'defaulttextinterpreter', 'tex')
%% Parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
lambda = [1 1.2];

%Con = [18 16];
Coff = [1 1];
M_int = [0 1; -1 1];
%K = [0 9; 11 6];

hill = Inf;
noise = 0;

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

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

% wave properties
num_waves = 1;
wave_type = 1; % because loading files may cause 
%% Load parameter sets for which travelling waves are predicted
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\TW_propagation_conditions_analytical';
fname_str = sprintf('trav_wave_conditions_check_wave_num_%d_type_%d_analysed_run_stability',...
    num_waves, wave_type);
load(fullfile(load_folder, fname_str), 'wave_possible', 'P');

% Display found waves
[x_found, y_found] = find(wave_possible);
t=table(P(x_found, :), y_found, 'VariableNames', {'Wave_type', 'Network'});
t2=table(x_found, y_found, 'VariableNames', {'wave_idx', 'Network'});
disp('Cell states F, M, B, E');
disp(t);
disp(t2);

% get list of all interaction matrices
M_int_found = get_found_M_int(y_found);
%% Run single simulation (test) for specific network/wave
% (Skip if testing for all types of waves)
%{
loop_idx = 1;

wave_idx = x_found(loop_idx);
network = y_found(loop_idx);
M_int = M_int_found{loop_idx};
states_perm = P(wave_idx, :); % new permutation
orig_perm = [2 4 3 1]; % default permutation F=(0,1)=2, M=(1,1)=4, B=(1,0)=3, E=(0,0)=1 => [2 4 3 1]

%==============================================================================
% Load predictor data
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\TW_propagation_conditions_analytical\run_stability';
fname_str_pred = sprintf('TW_conds_OR_logic_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
        num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));
load( fullfile(load_folder, fname_str_pred) );

%Con_wave = squeeze(Con_all(network, squeeze(trav_wave_conds_met(network, :)), :));
%K_wave = squeeze(K_all(network, squeeze(trav_wave_conds_met(network, :)), :, :));
Con_wave = squeeze( Con_all(network, :, :) );
K_wave = squeeze( K_all(network, :, :, :) );

%% ==============================================================================
% Load initial conditions with correct ordering
load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_snapshots';
wave_type_all_str = {'trav_wave_single_vertical',...
    'trav_wave_single_horizontal_inward_bend', 'trav_wave_single_horizontal_outward_bend'}; 
wave_type_str = wave_type_all_str{wave_type};

% load a default wave with default ordering
fname = fullfile(load_folder, wave_type_str);
cells_load = cell(2,1);
cells_load{1} = xlsread(fname, 'Sheet1');
cells_load{2} = xlsread(fname, 'Sheet2');
if all(size(cells_load{1})==[N 2])
    cells_in_temp = cells_load{1};
elseif all(size(cells_load{1})==[gz gz]) && all(size(cells_load{2})==[gz gz])
    cells_in_temp(:, 1) = reshape(cells_load{1}, N, 1);
    cells_in_temp(:, 2) = reshape(cells_load{2}, N, 1);
else
    disp('Wrong input format');
end

% get cell state population of initial state
cells_idx_temp = cells_in_temp*[2; 1];
%n_in = histcounts(cells_idx, -0.5:3.5);

% re-order cells according to states_perm for this wave
cells_in = zeros(N, 2);
%idx=cell(4,1);
states_default = {[0 0], [0 1], [1 0], [1 1]};
% adjust cells to have new permutation
for i=1:4
    idx_temp = find(cells_idx_temp==orig_perm(i)-1);
    perm = states_perm(i);
    cells_in(idx_temp, :) = repmat(states_default{perm}, numel(idx_temp), 1);    
end

hin=figure;
plot_handle = reset_cell_figure(hin, pos, rcell);
update_figure_periodic_scatter(plot_handle, cells_in, 0, 12, 0, a0, dist);
%% ==============================================================================
% Set parameters from predictor data
idx_waves = find(trav_wave_conds_met(network, :));
idx = idx_waves(27);
fprintf('Theoretical prediction: trav. wave possible? %d \n',...
    trav_wave_conds_met(network, idx) );
Con = Con_wave(idx,:);
K = squeeze(K_wave(idx, :, :));
M_int = M_int_found{loop_idx};
%%
% Default parameters
tmax = 100;
disp_mol = 12;
showI = 0;

% get cell state population of initial state
cells_idx = cells_in*[2; 1];
n_in = histcounts(cells_idx, -0.5:3.5);

% Run tester
display = 1;
[trav_wave, trav_wave_2, cellsOut] = trav_wave_sim_tester(cells_in, pos, dist, rcell, a0, M_int, Con,...
    Coff, K, lambda, hill, noise, tmax, display);
fprintf('Travelling wave? Method 1: %d (1=Yes, 0=No) \n', trav_wave);
fprintf('Travelling wave? Method 2: %d (1=Yes, 0=No) \n', trav_wave_2);

% check whether population has changed
cells_out_idx = cellsOut*[2; 1];
n_out = histcounts(cells_out_idx, -0.5:3.5);
fprintf('Unmodified? %d \n', all(n_out==n_in));
%}
%{
cells = cells_in;
cells_hist = {};
cells_hist{end+1} = cells;
    
hin=figure;
plot_handle = reset_cell_figure(hin, pos, rcell);
t = 0;
period = Inf; %default values
t_onset = Inf; 
[cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
        Rcell, Con, Coff, K, lambda, hill, noise);
update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);

% always check within first t_ac time steps
t_ac = 10^2; 
while changed && period==Inf && t<t_ac
    %disp(t);
    pause(0.1);
    t = t+1;
    cells = cellsOut;
    cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
    [period, t_onset] = periodicity_test_short(cells_hist); 
    update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);
    [cellsOut, changed] = update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
        Rcell, Con, Coff, K, lambda, hill, noise);
end

t_out = t; %default t_out
% if periodicity found, refine check to find period
if period<Inf && t>t_ac
    [period, t_onset] = periodicity_test_detailed(cells_hist, t_check, period);
    t_out = t_onset + period; 
end

if changed && t==tmax
    tmax_string = '_tmax_reached';
else
    tmax_string = '';
end
fprintf('Final: t_out = %d, period %d \n', t_out, period);

digits = 5;
[trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0, period, t_out, dist, digits);
fprintf('Travelling wave? Method 1: %d (1=Yes, 0=No) \n', trav_wave);
fprintf('Travelling wave? Method 2: %d (1=Yes, 0=No) \n', trav_wave_2);
%}
%% Check for a set of parameters
%
for loop_idx=1:numel(x_found)
    wave_idx = x_found(loop_idx);
    network = y_found(loop_idx);
    M_int = M_int_found{loop_idx};
    states_perm = P(wave_idx, :); % new permutation
    orig_perm = [2 4 3 1]; % default permutation F=(0,1)=2, M=(1,1)=4, B=(1,0)=3, E=(0,0)=1 => [2 4 3 1]
    fprintf('wave_idx = %d, network = %d \n', wave_idx, network);
    
    %==============================================================================
    % Load predictor data
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\TW_propagation_conditions_analytical\run_stability';
    fname_str_pred = sprintf('TW_conds_OR_logic_wave_num_%d_type_%d_states_%d_%d_%d_%d',...
        num_waves, wave_type, states_perm(1), states_perm(2), states_perm(3), states_perm(4));

    load( fullfile(load_folder, fname_str_pred) );
    Con_wave = squeeze(Con_all(network, :, :));
    K_wave = squeeze(K_all(network, :, :, :));
    fprintf('Sims to do = %d \n', size(Con_wave, 1));
    %}
    %% ==============================================================================
    % Load initial conditions with correct ordering
    load_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\travelling_wave_snapshots';
    wave_type_all_str = {'trav_wave_single_vertical',...
        'trav_wave_single_horizontal_inward_bend', 'trav_wave_single_horizontal_outward_bend'}; 
    wave_type_str = wave_type_all_str{wave_type};
    
    % load a default wave with default ordering
    fname = fullfile(load_folder, wave_type_str);
    cells_load = cell(2,1);
    cells_load{1} = xlsread(fname, 'Sheet1');
    cells_load{2} = xlsread(fname, 'Sheet2');
    if all(size(cells_load{1})==[N 2])
        cells_in_temp = cells_load{1};
    elseif all(size(cells_load{1})==[gz gz]) && all(size(cells_load{2})==[gz gz])
        cells_in_temp(:, 1) = reshape(cells_load{1}, N, 1);
        cells_in_temp(:, 2) = reshape(cells_load{2}, N, 1);
    else
        disp('Wrong input format');
    end
     
    %hin=figure;
    %plot_handle = reset_cell_figure(hin, pos, rcell);
    %update_figure_periodic_scatter(plot_handle, cells_in_temp, 0, 12, 0, a0, dist);
    %%
    % get cell state population of initial state
    cells_idx_temp = cells_in_temp*[2; 1];

    % re-order cells according to states_perm for this wave
    cells_in = zeros(N, 2);
    %idx=cell(4,1);
    states_default = {[0 0], [0 1], [1 0], [1 1]};
    % adjust cells to have new permutation
    for i=1:4
        idx_temp = find(cells_idx_temp==orig_perm(i)-1);
        perm = states_perm(i);
        cells_in(idx_temp, :) = repmat(states_default{perm}, numel(idx_temp), 1);    
    end

    % Display initial config
    %{
    qsave = 0;
    hin=figure;
    plot_handle = reset_cell_figure(hin, pos, rcell);
    update_figure_periodic_scatter(plot_handle, cells_in, 0, 12, 0, a0, dist);
    fname_str = sprintf('Wave_picture_network_%d_states_%d_%d_%d_%d', network,...
        states_perm(1), states_perm(2), states_perm(3), states_perm(4));
    fname_out = fullfile(save_folder, fname_str);
    save_figure(hin, 0, 0, fname_out, '.pdf', qsave, 1);
    
    continue
    %}
    %% ==============================================================================
    % Fixed parameters
    tmax = 100; % evolve system for tmax steps
    M_int = M_int_found{loop_idx};
    %disp_mol = 12;
    %showI = 0;
    display = 0;
    
    % variables to store
    n_pset_wave = size(K_wave, 1);
    trav_wave_all = zeros(n_pset_wave, 1); % is it still a trav wave?
    trav_wave_all_strict = zeros(n_pset_wave, 1); % strict method: p, I both constant (up to ... digits)
    unmodified_all = zeros(n_pset_wave, 1); % is the travelling wave unmodified (on population level)
    
    % get cell state population of initial state
    cells_idx = cells_in*[2; 1];
    n_in = histcounts(cells_idx, -0.5:3.5);

    for i=1:n_pset_wave
        disp(i);
        thisK = squeeze(K_wave(i, :, :));
        thisCon = Con_wave(i,:);
        %fprintf('K = [0 %d; %d %d] \n', thisK(1,2), thisK(2,1), thisK(2,2));
        
        [trav_wave, trav_wave_2, cells_out] = trav_wave_sim_tester(cells_in, pos, dist,...
            rcell, a0, M_int, thisCon, Coff, thisK, lambda, hill, noise, tmax, display);

        % check whether population has changed
        cells_out_idx = cells_out*[2; 1];
        n_out = histcounts(cells_out_idx, -0.5:3.5);
        fprintf('Unmodified? %d \n', all(n_out==n_in));

        % save results
        trav_wave_all(i) = trav_wave_2; % whether it is still a travelling wave
        trav_wave_all_strict(i) = trav_wave;
        unmodified_all(i) = all(n_out==n_in); % whether the population of the wave is unmodified
    end
    close all
    %}

    %% Save data
    %
    save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\TW_propagation_conditions_simulations';
    fname_str_data = sprintf('stability_sim_from_pred_all_wave_num_%d_type_%d_network_%d_states_%d_%d_%d_%d_OR_logic',...
         num_waves, wave_type, network, states_perm(1), states_perm(2), states_perm(3), states_perm(4));  
    fname = fullfile(save_folder, fname_str_data);
    save(fname, 'trav_wave_all', 'trav_wave_all_strict', 'unmodified_all',...
        'wave_idx', 'network', 'M_int', 'states_perm', 'fname_str_pred');

end
%}
%% functions
%[trav_wave] = trav_wave_sim_tester(cells_in, pos, dist, rcell, a0, M_int, Con,...
%    Coff, K, lambda, hill, noise, tmax);
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

function [trav_wave, trav_wave_2, cellsOut] = trav_wave_sim_tester(cells_in, pos, dist, rcell, a0, M_int, Con,...
    Coff, K, lambda, hill, noise, tmax, display)
    logic = 2; % OR logic
    disp_mol = 12;
    showI = 0;
    Rcell = rcell*a0;
    cells = cells_in;
    cells_hist = {};
    cells_hist{end+1} = cells;
    
    t = 0;
    period = Inf; %default values
    %t_onset = Inf; 
    [cellsOut, changed] = update_cells_two_signals_finite_Hill_new(cells, dist, M_int, a0,...
            Rcell, Con, Coff, K, lambda, hill, noise, logic);    
    if display
        hin=figure;
        plot_handle = reset_cell_figure(hin, pos, rcell);
        update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);
    end

    % always check within first t_ac time steps
    while changed && period==Inf && t<tmax
        %disp(t);
        pause(0.1);
        t = t+1;
        cells = cellsOut;
        cells_hist{end+1} = cells; %{cells(:, 1), cells(:, 2)};
        [period, ~] = periodicity_test_short(cells_hist); 
        if display
            pause(0.2);
            update_figure_periodic_scatter(plot_handle, cells, t, disp_mol, showI, a0, dist);
        end
        [cellsOut, changed] = update_cells_two_signals_finite_Hill_new(cells, dist, M_int, a0,...
        	Rcell, Con, Coff, K, lambda, hill, noise, logic);    
    end

    t_out = t;
    fprintf('Final: t_out = %d, period %d \n', t_out, period);
    
    % check whether produced pattern is still a travelling wave
    digits = 5;
    [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0, period, t_out, dist, digits);

    fprintf('Travelling wave? Meth. 1: %d (1=Yes, 0=No) \n', trav_wave);
    fprintf('Travelling wave? Meth. 2: %d (1=Yes, 0=No) \n', trav_wave_2);
end