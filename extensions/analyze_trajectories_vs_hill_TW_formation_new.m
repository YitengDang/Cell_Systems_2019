%% Analyze saved trajectories across a range of parameters
% new: new classification based on homogeneity and stationarity
clear all
close all
set(0, 'defaulttextinterpreter', 'tex');

%% Parameters
N = 225;
tmax = 10000;
hill_all = [Inf 100 50 20 10 5:-1:1];
network = 15;

var_all = hill_all; %mcsteps_all;
label1 = 'vs_hill';
label2 = 'network_15';

num_params = 10;
nruns = 100; %number of runs per parameter set
vpa_digits = 36;

% folder for saving figures
save_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_hill';
subfolder = sprintf('TW_formation_network_%d', network);
save_path_fig = fullfile(save_path, subfolder);
%% Get TW frequencies from loaded data
filecount = zeros(numel(var_all), num_params); 
t_out_full = zeros(numel(var_all), num_params, nruns); % final times
period_full = zeros(numel(var_all), num_params, nruns); % periodicity test
t_onset_full = zeros(numel(var_all), num_params, nruns); 
trav_wave_full = zeros(numel(var_all), num_params, nruns); 
trav_wave_full_2 = zeros(numel(var_all), num_params, nruns); 

%% Load negative control (no noise, infinite Hill)
load_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_reliability';
load_file = 'analyzed_data_TW_formation_network_15_nruns_500_digits_5';
%load_file = 'analyzed_data_TW_formation_network_15_nruns_500_digits_5_with_K_Con_data';
temp = load(fullfile(load_path, load_file));

% store loaded data
t_out_full(1, :, :) = temp.t_out_all(1:num_params, 1:nruns);
period_full(1, :, :) = temp.period_all(1:num_params, 1:nruns);
t_onset_full(1, :, :) = temp.t_onset_all(1:num_params, 1:nruns);
trav_wave_full(1, :, :) = temp.trav_wave_all(1:num_params, 1:nruns);
trav_wave_full_2(1, :, :) = temp.trav_wave_all_2(1:num_params, 1:nruns);
%% Load data files
%
%%% PART I %%%
% folders for loading data
load_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_hill';
subfolder = sprintf('TW_formation_network_%d', network);

% store parameter sets
Con_all_sim = zeros(num_params, 2); 
K_all_sim = zeros(num_params, 2, 2);

for p_set = 1:num_params
    disp(p_set);
    subsubfolder = sprintf('Parameter_set_%d', p_set);
    folder = fullfile(load_path, subfolder, subsubfolder);
    
    % --- Manually set ---
    pattern = sprintf(...
        'two_signal_mult_N%d_params_%s_sim_%s_hill_%s_vpa_digits_%d_t_out_%s_period_%s_vpa_digits_%d-v%s',...
        N, '(\d+)', '(\d+)', '(\d+)', vpa_digits, '(\d+)', '(\d+|Inf)',...
        vpa_digits, '(\d+)'); % '.' = anything
    
    % get all parameters from first simulation
    listing = dir(folder);
    num_files = numel(listing)-2; %first two entries are not useful
    for i = 1:num_files
        filename = listing(i+2).name;
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat') && ~isempty(regexp(filename, pattern, 'once'))
            load( fullfile( folder, filename), 'save_consts_struct');
            Con_all_sim(p_set, :) = save_consts_struct.Con;
            K_all_sim(p_set, :, :) = save_consts_struct.K;
            break
        end
    end
    
    % get filenames of all simulations
    names = {};
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
    
    % -------------------- 
    % load all simulations    
    for i=1:numel(names)
        if isempty(regexp(names{i}, pattern, 'once')) % only load files matching a certain pattern
            continue
        else
            disp(names{i});
            load( fullfile( folder, strcat(names{i}, '.mat')), 'cells_hist',...
                'save_consts_struct', 'distances', 'period', 't_out', 't_onset');

            [tokens, ~] = regexp(names{i}, pattern, 'tokens', 'match');

            % --- Manually set ---
            hill = save_consts_struct.hill;
            idx = find(hill == hill_all, 1);
            % --------------------
            %{
            if num_params>1
                idx2 = str2double(tokens{1}{1});
            elseif num_params==1
                idx2 = 1;
            end        
            %}
            idx2 = p_set;
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
                
                % store results
                t_out_full(idx, idx2, idx3) = t_out;
                period_full(idx, idx2, idx3) = period;
                t_onset_full(idx, idx2, idx3) = t_onset;

                % test for TW
                a0 = save_consts_struct.a0;
                digits = 3;
                [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                    period, t_out, distances, digits);
                trav_wave_full(idx, idx2, idx3) = trav_wave;
                trav_wave_full_2(idx, idx2, idx3) = trav_wave_2;
            end
            %
        end
    end
    %}
end

%{
t_out_full = t_out_full(1:5,:,:);
period_full = period_full(1:5,:,:);
t_onset_full = t_onset_full(1:5,:,:);
trav_wave_full = trav_wave_full(1:5,:,:);
trav_wave_full_2 = trav_wave_full_2(1:5,:,:);
%}
%% Save the loaded data
%
fname_str = sprintf('analyzed_data_%s_vs_hill_num_params_%d_nruns_%d_digits_%d',...
    subfolder, num_params, nruns, digits);
save_data_path = fullfile(save_path, subfolder);
%save( fullfile(save_data_path, strcat(fname_str, '.mat')), 'noise_all',...
%    'filecount', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
%    'num_params', 'nruns', 'digits', 'trav_wave_all', 'trav_wave_all_2');
save( fullfile(save_data_path, strcat(fname_str, '.mat')), 'hill_all',...
    'filecount', 't_out_full', 'period_full', 't_onset_full', 'tmax',...
    'num_params', 'nruns', 'vpa_digits', 'trav_wave_full', 'trav_wave_full_2');
%% Load the saved data
fname_str = sprintf('analyzed_data_%s_vs_hill_num_params_%d_nruns_%d_digits_%d',...
    subfolder, num_params, nruns, digits);
data_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_hill\TW_formation_network_15';
%save( fullfile(save_data_path, strcat(fname_str, '.mat')), 'noise_all',...
%    'filecount', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
%    'num_params', 'nruns', 'digits', 'trav_wave_all', 'trav_wave_all_2');
load( fullfile(data_path, strcat(fname_str, '.mat')), 'hill_all',...
    'filecount', 't_out_full', 'period_full', 't_onset_full', 'tmax',...
    'num_params', 'nruns', 'vpa_digits', 'trav_wave_full', 'trav_wave_full_2');

%% 4-way classification
% (1) static, (2) oscillatory, (3) TW end states, (4) infinite dynamics
%
idx3 = (t_out_full < tmax & period_full == Inf);
idx2 = (t_out_full < tmax & period_full < Inf & ~trav_wave_full_2);
idx1 = (trav_wave_full_2);
idx4 = (t_out_full == tmax);
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

%% Plot fractions according to 4-way classification
h = figure;
x_data = var_all;

% x-axis custom scale
bar(flipud(frac_all), 'stacked');
xticks = fliplr(hill_all); xticks(end) = 10^3;
xtick_labels = sprintfc('%d', fliplr(hill_all(1:end)));
set(gca, 'XTick', 1:numel(hill_all), 'XTickLabels', xtick_labels );

% x-axis log-scale 
%bar( log10(x_data), frac_all, 'stacked');
%set(gca, 'XTick', -3:0, 'XTickLabels', sprintfc('10^{%d}', -3:0) );

% x-axis evenly spread
%bar( frac_all, 'stacked');
%set(gca, 'XTick', 1:numel(x_data), 'XTickLabels', sprintfc('%.3f', x_data) );

xlabel('Hill coefficient');
%xlabel('Noise strength \alpha/K^{(ij)}');
%xlabel('Lattice disorder (Monte Carlo steps)');
ylabel('Fraction of simulations');
set(gca, 'FontSize', 32);
box on
% classificiation I
legend({'Travelling wave', 'Oscillatory', 'Static', 'Infinite dynamics'}, 'location', 'sw');
% classification Douwe
%legend({'Static homogeneous', 'Homogeneous oscillations', 'Pure wave', 'Infinite dynamics'});
%legend({'Pure travelling wave', 'Infinite dynamics', 'Homogeneous, oscillatory', 'Homogeneous, static'});

qsave = 1;
fname = fullfile(save_path_fig, sprintf('TW_formation_%s_%s_nruns_%d_digits_%d_%s',...
    label1, label2, nruns, digits, 'classification_dynamics'));
save_figure(h, 12, 8, fname, '.pdf', qsave);

%% Analyze fraction of TWs
%trav_wave_all_mean = sum(sum(trav_wave_all, 3), 2)/(num_params*nruns);
trav_wave_all_2_mean = sum(sum(trav_wave_full_2, 3), 2)/(num_params*nruns);

h = figure;
x_data = hill_all;    
x_data(1) = 10^(2.5);
hold on
plot( x_data, trav_wave_all_2_mean, 'k--', 'LineWidth', 1.5 );
scatter(x_data, trav_wave_all_2_mean, 'k', 'filled')
ylim([0 0.3]);
set(gca, 'XScale', 'log');
box on

% labels (change per case)
xlabel('Hill coefficient');
%xlabel('Lattice disorder (Monte Carlo steps)');
%xlabel('Noise strength \alpha/K^{(ij)}');
ylabel('Fraction TW');

% Ticks
xlim([1 10^2.5]);
xticks = 10.^([0:2 2.5]);
% Tick labels
xtick_labels = sprintfc('10^{%d}', 0:3);
xtick_labels{end} = 'Inf';
set(gca, 'FontSize', 32)
set(gca, 'XTick', xticks, 'XTickLabels', xtick_labels);
set(gca, 'xdir', 'reverse') % reverse axes direction

qsave = 1;
fname = fullfile(save_path_fig, sprintf('TW_formation_%s_%s_nruns_%d_digits_%d_%s',...
    label1, label2, nruns, digits, 'frac_TW_all_10_8_final'));
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Analyze fraction of TWs per parameter set
trav_wave_all_2_mean = sum(trav_wave_full_2, 3)/(nruns);

h = figure;
x_data = hill_all;
x_data(1) = 10^3;
plot( x_data, trav_wave_all_2_mean, 'o-', 'LineWidth', 1.5 );
ylim([0 1]);
set(gca, 'XScale', 'log');

% labels (change per case)
xlabel('Hill coefficient');
%xlabel('Lattice disorder (Monte Carlo steps)');
%xlabel('Noise strength \alpha/K^{(ij)}');
ylabel('Fraction TW');

% Ticks
xticks = 10.^(0:3);
% Tick labels
xtick_labels = sprintfc('10^{%d}', 0:3);
xtick_labels{end} = 'Inf';
set(gca, 'FontSize', 32, 'XTick', xticks, 'XTickLabels', xtick_labels);
set(gca, 'xdir', 'reverse') % reverse axes direction

% legend
legend(sprintfc('%d', 1:10))
qsave = 1;
fname = fullfile(save_path_fig, sprintf('TW_formation_%s_%s_nruns_%d_digits_%d_%s',...
    label1, label2, nruns, digits, 'frac_TW_all_by_pset_10_8'));
save_figure(h, 12, 8, fname, '.pdf', qsave);
%% Analyze max t_out
% # trajectories reaching tmax
t_max_reached = (t_out_full == tmax);
t_max_count = sum(sum(t_max_reached, 3), 2);

h = figure;
x_data = hill_all;
x_data(1) = 10^3;
plot( x_data, t_max_count/(num_params*nruns), 'bo-', 'LineWidth', 1.5  );
%plot( x_data, t_out_mean, 'bo-', 'LineWidth', 1.5  );
%errorbar( x_data, t_out_mean, t_out_var, 'bo-', 'LineWidth', 1.5  );
%xlim([x_data(1) x_data(end)]);
ylim([0 1]);
set(gca, 'XScale', 'log');
xlabel('Hill coefficient');
ylabel('Fraction reaching t_{max}');
title(sprintf('t_{max}=10^{%d}', log10(tmax)));

%imagesc(p1, p2, t_max_count/nruns);
%set(gca, 'YDir', 'Normal');
%c = colorbar;
%caxis([0 1]);
%colormap('parula');
%xlabel('$$p_1$$')
%ylabel('$$p_2$$')
%ylabel(c, 'fraction');

% Ticks
xticks = 10.^(0:3);
% Tick labels
xtick_labels = sprintfc('10^{%d}', 0:3);
xtick_labels{end} = 'Inf';
set(gca, 'FontSize', 32, 'XTick', xticks, 'XTickLabels', xtick_labels);
set(gca, 'xdir', 'reverse') % reverse axes direction

qsave = 1;
fname = fullfile(save_path_fig, sprintf('TW_formation_%s_%s_nruns_%d_digits_%d_%s',...
    label1, label2, nruns, digits, 'frac_reaching_tmax'));
save_figure(h, 12, 8, fname, '.pdf', qsave);

%% Analyze mean t_out
% (1) <t_out>
t_out_mean = mean(mean(t_out_full, 3), 2);
t_out_std = std(t_out_full(:,:), 0, 2);

h = figure;
x_data = hill_all;
x_data(1) = 10^3;
%plot( x_data, t_out_mean, 'bo-', 'LineWidth', 1.5  );
errorbar( x_data, t_out_mean, t_out_std, 'bo-', 'LineWidth', 1.5  );
%xlim([x_data(1) x_data(end)]);
%ylim([0 1]);
set(gca, 'XScale', 'log');
xlabel('Hill coefficient');
ylabel('\langle t_{out} \rangle');
title(sprintf('t_{max}=10^{%d}', log10(tmax)));

%imagesc(p1, p2, t_max_count/nruns);
%set(gca, 'YDir', 'Normal');
%c = colorbar;
%caxis([0 1]);
%colormap('parula');
%xlabel('$$p_1$$')
%ylabel('$$p_2$$')
%ylabel(c, 'fraction');

% Ticks
xticks = 10.^(0:3);
% Tick labels
xtick_labels = sprintfc('10^{%d}', 0:3);
xtick_labels{end} = 'Inf';
set(gca, 'FontSize', 32, 'XTick', xticks, 'XTickLabels', xtick_labels);
set(gca, 'xdir', 'reverse') % reverse axes direction

qsave = 1;
fname = fullfile(save_path_fig, sprintf('TW_formation_%s_%s_nruns_%d_digits_%d_%s',...
    label1, label2, nruns, digits, 't_out_mean'));
save_figure(h, 12, 8, fname, '.pdf', qsave);
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