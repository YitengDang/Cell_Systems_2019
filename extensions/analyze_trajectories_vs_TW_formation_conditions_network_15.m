%% Analyze saved trajectories across a range of parameters
clear all
close all
set(0, 'defaulttextinterpreter', 'tex');

%% Parameters
N = 225;
tmax = 10000;
network = 15;

num_params = 2534;
nruns = 10; %number of runs per parameter set

remote = 0;

% folder for saving figures
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_formation_with_prop_conditions';
%% Load data files
%{
%%% PART I %%%
% folders for loading data
load_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\TW_formation_with_prop_parameters';

subfolder = sprintf('TW_formation_network_%d', network);
folder = fullfile(load_path, subfolder);

listing = dir(folder);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
    end
end

%% Load data
filecount = zeros(1, num_params); 
t_out_all = zeros(1, num_params, nruns); % final times
period_all = zeros(1, num_params, nruns); % periodicity test
t_onset_all = zeros(1, num_params, nruns); 
trav_wave_all = zeros(1, num_params, nruns); 
trav_wave_all_2 = zeros(1, num_params, nruns); 
hom_end_state_all = zeros(1, num_params, nruns); % whether end state is homogeneous

% --- Manually set ---
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_rand_params_%s_t_out_%s_period_%s-v%s',...
    N, '(\d+)', '(\d+)', '(\d+|Inf)', '(\d+)'); % '.' = anything
% --------------------

for i=1:numel(names)
    if isempty(regexp(names{i}, pattern, 'once')) % only load files matching a certain pattern
        continue
    else
        %disp(names{i});
        load( fullfile( folder, strcat(names{i}, '.mat')), 'cells_hist',...
            'save_consts_struct', 'distances', 'period', 't_out', 't_onset');
        
        [tokens, ~] = regexp(names{i}, pattern, 'tokens', 'match');
        
        % --- Manually set ---
        %mcsteps = save_consts_struct.mcsteps;
        %idx = find(mcsteps == mcsteps_all, 1);
        
        %noise = save_consts_struct.noise;
        idx = 1; %find(noise == noise_all, 1);
        % --------------------
        
        idx2 = str2double(tokens{1}{1});
        %disp(idx2);
        %
        if ~isempty(idx)
            filecount(idx, idx2) = filecount(idx, idx2) + 1;
            idx3 = filecount(idx, idx2);
            if idx3>nruns
                % only do up to nruns simulations
                continue
            end
            %disp(idx3);
            disp(names{i});
            
            t_out_all(idx, idx2, idx3) = t_out;
            period_all(idx, idx2, idx3) = period;
            t_onset_all(idx, idx2, idx3) = t_onset;
            
            % test for TW
            a0 = save_consts_struct.a0;
            digits = 3;
            [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                period, t_out, distances, digits);
            trav_wave_all(idx, idx2, idx3) = trav_wave;
            trav_wave_all_2(idx, idx2, idx3) = trav_wave_2;
            
            % check end state
            cells_final = cells_hist{end};
            hom_end_state = size(unique(cells_final, 'rows'), 1)==1;
            hom_end_state_all(idx, idx2, idx3) = hom_end_state;
        end
        %
    end
end
%}
%% Save the loaded data
%{
subfolder = sprintf('TW_propagation_network_%d', network);
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);

save_data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\TW_formation_with_prop_parameters';
save( fullfile(save_data_path, strcat(fname_str, '.mat')),...
    'filecount', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
    'num_params', 'nruns', 'digits', 'trav_wave_all', 'trav_wave_all_2',...
    'hom_end_state_all');
%}
%% Load the saved data
%
subfolder = sprintf('TW_propagation_network_%d', network);
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);

load_data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\TW_formation_with_prop_parameters';
load( fullfile(load_data_path, strcat(fname_str, '.mat')),...
    'filecount', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
    'num_params', 'nruns', 'digits', 'trav_wave_all', 'trav_wave_all_2',...
    'hom_end_state_all');
%}
%% 4-way classification
%
% (1) static, (2) oscillatory, (3) TW end states, (4) infinite dynamics
%
idx3 = (t_out_all < tmax & period_all == Inf);
idx2 = (t_out_all < tmax & period_all < Inf & ~trav_wave_all_2);
idx1 = (trav_wave_all_2);
idx4 = (t_out_all == tmax);
% all(all(all(idx1 + idx2 + idx3 + idx4)))
%}

% Douwe's classification (4) static homogeneous, (3) oscillatory homogeneous, (1) TW end states, (2) infinite dynamics
%{
idx1 = (trav_wave_all_2);
idx3 = (t_out_all < tmax & period_all == Inf & hom_end_state_all);
idx4 = (t_out_all < tmax & period_all < Inf & ~trav_wave_all_2 & hom_end_state_all);
idx2 = ones(size(idx1)) - (idx1+idx4+idx3); %(t_out_all == tmax);
%all(all(all(idx1 + idx2 + idx3 + idx4)))
%}
% calculate fractions
frac_all = [sum(sum(idx1, 3), 2)/(num_params*nruns) sum(sum(idx2, 3), 2)/(num_params*nruns)...
    sum(sum(idx3, 3), 2)/(num_params*nruns) sum(sum(idx4, 3), 2)/(num_params*nruns)];
frac_all_by_pset = [sum(idx1, 3)/(nruns); sum(idx2, 3)/(nruns);...
    sum(idx3, 3)/(nruns); sum(idx4, 3)/(nruns)];

%% TW formation per parameter set
% fraction of parameter sets generating at least one wave
num_wave_pset = sum(frac_all_by_pset(1,:)>0);
frac_wave_pset = num_wave_pset/num_params;

% distribution of fractions of TWs per parameter set (mean + std)
h = figure;
histogram(frac_all_by_pset(1,:), -0.05:0.1:1.05)
%histogram(frac_all_by_pset(1,frac_all_by_pset(1,:)>0), -0.05:0.1:1.05)
title(sprintf('%d/%d TW formation', num_wave_pset, num_params));
xlabel('Fraction of simulations forming TW');
ylabel('Number of parameter sets');
set(gca, 'FontSize', 32);
box on

qsave = 1;
fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
    sprintf('_nruns_%d_digits_%d', nruns, digits), '_fraction_TW_formed_by_pset_hist'));
save_figure(h, 12, 8, fname, '.pdf', qsave);
 

%% Analyze further which parameter sets give higher TW formation chances
TW_param_idx = frac_all_by_pset(1,:)>0; % indices of parameter sets giving rise to waves

% Load simulation parameters
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2_net_parameters_TW_sim';
if remote
    folder = strrep(folder, 'N:\', 'W:\staff-bulk\');
end
fname_str = 'Wave_type_1_network_15_states_F3_M4_B2_E1_Con_K_values_waves_sim';
load(fullfile(folder, fname_str), 'N', 'a0', 'hill', 'lambda', 'noise',...
    'rcell', 'Con_wave_sim', 'K_wave_sim');
M_int = [1 -1; 1 0];

%% Plot distributions of parameters

% Plot spider plots
% filter data based on present interactions (not generalized)
K_idx = 1:3;

% input vars
P_data = [Con_wave_sim(TW_param_idx, :) K_wave_sim(TW_param_idx, K_idx)];
P_labels = {'C_{ON}^{(1)}', 'C_{ON}^{(2)}', 'K^{(11)}',...
        'K^{(12)}', 'K^{(21)}', 'K^{(22)}'};
axes_interval = 3;

% plot
range_values_in = repmat([1 1000], 5, 1);
spider_plot(P_data, P_labels([1:2 K_idx+2]), axes_interval, range_values_in,...
    'Marker', 'o',...
    'LineStyle', '-',...
    'Color', [1 0 0 0.02],...
    'LineWidth', 2,...
    'MarkerSize', 2);

%spider_plot_linear(P_data, P_labels([1:2 K_idx+2]), axes_interval);
%title(sprintf('n=%d, nw %d', size(P_data,1), network),...
%    'Fontweight', 'bold', 'FontSize', 28);
set(gcf, 'Units', 'Inches', 'Position', [2 2 10 8]);
h = gcf;
hold off 

% save plot
set(h, 'Units', 'Inches', 'Position', [0.1 0.1 10 8]);
qsave = 1;
fname_str = sprintf('Parameters_TW_spider_plot_network_%d_linear_TW_formed', network);
save_figure(h, 10, 8, fullfile(save_path_fig, fname_str),'.pdf', qsave);

%% Plot parameters for which waves were not formed
h = figure;

% input vars
P_data = [Con_wave_sim(~TW_param_idx, :) K_wave_sim(~TW_param_idx, K_idx)];
P_labels = {'C_{ON}^{(1)}', 'C_{ON}^{(2)}', 'K^{(11)}',...
        'K^{(12)}', 'K^{(21)}', 'K^{(22)}'};
spider_plot(P_data, P_labels([1:2 K_idx+2]), axes_interval, range_values_in,...
    'Marker', 'o',...
    'LineStyle', '-',...
    'Color', [1 0 0 0.02],...
    'LineWidth', 2,...
    'MarkerSize', 2);
set(gcf, 'Units', 'Inches', 'Position', [2 2 10 8]);
hold off 

% save plot
set(h, 'Units', 'Inches', 'Position', [0.1 0.1 10 8]);
qsave = 1;
fname_str = sprintf('Parameters_TW_spider_plot_network_%d_linear_TW_not_formed', network);
save_figure(h, 10, 8, fullfile(save_path_fig, fname_str),'.pdf', qsave);

% histograms
%{
edges = 0:50:1000;

figure;
hold on
histogram( Con_wave_sim(TW_param_idx, 1), edges, 'normalization', 'pdf' );
histogram( Con_wave_sim(~TW_param_idx, 1), edges, 'normalization', 'pdf' );

figure;
hold on
histogram( Con_wave_sim(TW_param_idx, 2), edges, 'normalization', 'pdf' );
histogram( Con_wave_sim(~TW_param_idx, 2), edges, 'normalization', 'pdf' );

figure;
hold on
histogram( K_wave_sim(TW_param_idx, 1), edges, 'normalization', 'pdf' );
histogram( K_wave_sim(~TW_param_idx, 1), edges, 'normalization', 'pdf' );

figure;
hold on
histogram( K_wave_sim(TW_param_idx, 2), edges, 'normalization', 'pdf' );
histogram( K_wave_sim(~TW_param_idx, 2), edges, 'normalization', 'pdf' );

figure;
hold on
histogram( K_wave_sim(TW_param_idx, 3), edges, 'normalization', 'pdf' );
histogram( K_wave_sim(~TW_param_idx, 3), edges, 'normalization', 'pdf' );
%}

%% TW detection algorithm test
%{
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\randomized lattice\TW_formation_network_15';
fname_str = 'two_signal_mult_N225_ini_state_rand_params_29_mcsteps_0_t_out_2707_period_225-v1';

fname = fullfile(folder, fname_str);
load(fname, 'cells_hist', 'save_consts_struct', 'distances', 'period', 't_out');

a0 = save_consts_struct.a0;
digits = 3;
[trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
    period, t_out, distances, digits);
%}
%%
% (3) <t_out> over trajectories that do not reach t_max
%{
t_out_mean_2 = zeros( numel(mcsteps_all), 1 );
for i=1:numel(t_out_all)
    if ~t_max_reached(i)
        [i1,i2,i3] = ind2sub(size(t_out_all), i); 
        t_out_mean_2(i1, i2) = t_out_mean_2(i1, i2) + t_out_all(i1,i2,i3);
    end
end
t_out_mean = t_out_mean_2./(nruns - t_max_count);

h3 = figure(3);

set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, '$$\langle t_{out} |  t_{out} < t_{max} \rangle$$', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);
caxis([0 tmax]);
%ylim(c, [0 tmax]);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_mean_t_out_subset'));
    save_figure(h3, 10, 8, fname, '.pdf');
end
%}

%% Analyze periods
%{
% number of chaotic trajectories
n_chaotic = sum((period_all(:)==Inf).*(t_out_all(:) == tmax));

% data
period_data = period_all(period_all~=0);
uniq = unique(period_data);
C = categorical(period_data, uniq);

n_periodic = sum(period_data(:)~=Inf);

% distribution of periods
h4 = figure(4);
histogram(C);
xlabel('Period');
ylabel('Count');
title(sprintf('%d trajectories, %d periodic, %d chaotic',...
    sum(sum(filecount)), n_periodic, n_chaotic));
set(gca, 'FontSize', 20);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_periods_hist_chaotic'));
    save_figure(h4, 10, 8, fname, '.pdf');
end

%% average period (count only found periods)
period_mean = zeros(numel(p1), numel(p2));
period_count = zeros(numel(p1), numel(p2));
for i=1:numel(t_out_all)
    if period_all(i)>0 && period_all(i)<Inf
        [i1,i2,i3] = ind2sub(size(period_all), i); 
        period_count(i1, i2) = period_count(i1, i2) + 1;
        period_mean(i1, i2) = period_mean(i1, i2) + period_all(i1,i2,i3);
    end
end
period_mean = period_mean./period_count;

h5 = figure(5);
imagesc(p1, p2, period_mean);
set(gca, 'YDir', 'Normal');
c = colorbar;
caxis([0 max(uniq(uniq<Inf))])
xlabel('$$p_1$$')
ylabel('$$p_2$$')
ylabel(c, 'Mean period (time steps)', ...
    'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 20);

qsave = 0;
if qsave
    fname = fullfile(save_path_fig, strcat(fname_str, '_mean_period'));
    save_figure(h5, 10, 8, fname, '.pdf');
end
%}