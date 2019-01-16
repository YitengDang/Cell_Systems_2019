%% Plots one of the functions ("wave energy") devised to study "waveness" of a simulation over time
clear all
close all
clc
set(0,'defaulttextinterpreter', 'tex');

%% Import data
%folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\filmstrip_selection';
%fname_str = 'travelling_wave_osc_background_network_20_broad_TW_vertical';
folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\Fig5-self-organisation';
load_fname_str = 'sample_complex_trajectory_network_34';

load(fullfile(folder, load_fname_str), 'cells_hist', 'positions', 'distances', 'save_consts_struct');
rcell = save_consts_struct.rcell;
a0 = save_consts_struct.a0;
t_out = numel(cells_hist)-1;

% SPECIFY NETWORK
network = 34;
appendix = ''; % Note special rule for 33a, 33b
networks_all = [15 19 33 33 34 36];
network_idx = find(network==networks_all, 1);
if strcmp(appendix, 'b')
    network_idx = 4; % special case: 33b
end

[period_ub, t_onset] = periodicity_test_short(cells_hist);
[period, t_onset] = periodicity_test_detailed(cells_hist, t_onset,...
    period_ub);

%% get wave states
cell_states_all = cell(6, 1);
cell_states_all{1} = [1 0; 1 1; 0 1; 0 0]; % F, M, B, E
cell_states_all{2} = [1 1; 1 0; 0 0; 0 1]; % F, M, B, E
cell_states_all{3} = cell_states_all{1}; % F, M, B, E
cell_states_all{4} = [1 1; 0 1; 0 0; 1 0]; % F, M, B, E
cell_states_all{5} = cell_states_all{4}; % F, M, B, E
cell_states_all{6} = [0 1; 1 1; 1 0; 0 0]; % F, M, B, E
cell_states = cell_states_all{network_idx};

%% Compute wave energy (single state)
% wave energy parameters
Jself = 1;


% compute for given system state
cells = cells_hist{1};

