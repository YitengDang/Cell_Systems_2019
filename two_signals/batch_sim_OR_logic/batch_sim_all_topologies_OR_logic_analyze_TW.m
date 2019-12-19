%% Analyze data from batch simulations with OR logic
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'tex');

%% Load analyzed data
save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\batch_sim_all_topologies_OR_logic_run1\analyzed_data';
%save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\batch_sim_all_topologies_OR_logic_run1\analyzed_data_1sim';
fname_str = 'batch_sim_OR_logic_all_topologies_analyzed_all';
%fname_str = 'batch_sim_OR_logic_all_topologies_analyzed_network_3';
fname_out = fullfile(save_folder, fname_str);
load(fname_out);

%% Basic parameters
n_pset = 10^4;
nsim = 10;
tmax = 10^4; 
save_folder_fig = fullfile('H:\My Documents\Multicellular automaton\figures', ...
    'batch_sim_all_topologies_OR_logic_run1');

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
digits = -floor(log10(1/N)); % number of significant digits required to determine constancy of p

% get pos, dist
mcsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Network classification
% Based on following results, propose classification
class = cell(3, 1);
class{1} = [15	16	19	20	33	34	36	43]; % complex dynamics
class{2} = [2	5	8	10	11	12	13	14	18	21	22	23	25	27, ...
                29	30	31	32  35	37	38	39	40	41	42	44]; % simple oscillations
class{3} = [1	3	4	6	7	9	17	24	26	28]; % no oscillations

%% Plot fraction of TWs per network
frac_TW = sum(TW_count_loose_network_all, 2)/n_pset/nsim;
h = figure;
bar(1:numel(network_all), frac_TW);
xlabel('Network');
ylabel('Fraction TWs');
set(gca, 'FontSize', 32);
box on

% Save figures
set(h, 'Units', 'Inches', 'Position', [0.1 0.1 10 8]);
qsave = 0;
fname_str = 'Robustness_frac_sim_TW_all';
save_figure(h, 10, 8, fullfile(save_folder_fig, fname_str),'.pdf', qsave);

%% Plot unnormalized results
networks_idx = [15 19 33 34 36]; %class{1}; 
num_params = [5 5 6 6 6]; % number of parameters per network

% get # parameter sets that gave at least one TW
y_data = sum(TW_count_loose_network_all(networks_idx, :), 2)/n_pset/nsim;

% overall fraction of simulations that gave TWs
%y_data = TW_count_loose(networks_idx2)/(n_pset*nsim);

h = figure;
bar( 1:numel(networks_idx), y_data);
%set(gca, 'XTickLabel', sprintfc('Network %d', networks_sel(networks_idx2)) );
xlabel('Network');
set(gca, 'XTickLabel', sprintfc('%d', networks_idx) );
xtickangle(45)
%ylabel('Count');
ylabel('Frequency');
%set(gca, 'YScale', 'log');
set(gca, 'FontSize', 32);
%ylim([10^(-4) 1]);
%ylim([0 2*10^(-3)]);
box on

% Save figures
set(h, 'Units', 'Inches', 'Position', [0.1 0.1 10 8]);
qsave = 0;
fname_str = 'Robustness_frac_sim_TW_selected';
save_figure(h, 10, 8, fullfile(save_folder_fig, fname_str),'.pdf', qsave);

%% Filter out TW simulations from larger set
% Move files to a new folder
%{
save_folder_TW = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\batch_sim_all_topologies_OR_logic_run1\selected_sims\trav_waves';
save_folder_data = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals_batch_sim_2\batch_sim_all_topologies_OR_logic_run1';
[idx1, idx2] = find( TW_count_loose_network_all>0 );

for i=1:numel(idx1)
    %%
    network_old = network_all(idx1(i));
    network_new = idx1(i);
    pset = idx2(i);
    fprintf('Network %d, pset %d \n', network_new, pset);
    
    folder = fullfile(save_folder_data, sprintf('Network_%d', network_old), ...
        sprintf('Param_%d', pset) );
    all_files = dir(folder);
    for i2 = 3:numel(all_files)-1
        fname = all_files(i2).name;
        load(fullfile(folder, fname), 'cells_hist', 'save_consts_struct',...
            'period', 't_out', 'distances');
        gz = sqrt(save_consts_struct.N);
        a0 = save_consts_struct.a0;
        
        [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
            period, t_out, distances, digits);
        
        if trav_wave_2 
            % copy file 
            new_file_name = sprintf('Network_%d_Param_%d_%s', network_new, pset, fname);
            new_file = fullfile(save_folder_TW,  new_file_name);
            disp(new_file);
            
            old_file = fullfile(folder, fname);
            copyfile(old_file, new_file)
        end
    end
end
%}
%% Functions
function M_out = M_perm(M, perm)
    if perm==1
        M_out = M;
    elseif perm==2
        M_out = [M(2,2) M(2,1); M(1,2) M(1,1)];
    else
        warning('Wrong permutation');
        M_out = Inf;
    end
end