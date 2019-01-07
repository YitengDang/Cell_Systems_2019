%% Analyze saved trajectories across a range of parameters
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');

%% Parameters
N = 225;
tmax = 10000;
%mcsteps_all = [0 10 100 1000];
%mcsteps_all = [10 20:20:100 200:200:1000];
noise_all = [0.01 0.03 0.05 0.1 0.3 0.5 1];
noise_all_2 = [3 5 10 50 100];
network = 15;

var_all = noise_all;

num_params = 100;
nruns = 20; %number of runs per parameter set

% folder for saving figures
%save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_vs_mcsteps';
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_vs_noise';

%% Load data files 
%{
% (Part 1)
%path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\randomized lattice';
path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_with_noise';
subfolder = sprintf('TW_propagation_network_%d', network);

folder = fullfile(path, subfolder);

listing = dir(folder);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
names = {};
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
% (Part 1)
filecount = zeros(numel(var_all), num_params); 
t_out_all = zeros(numel(var_all), num_params, nruns); % final times
period_all = zeros(numel(var_all), num_params, nruns); % periodicity test
t_onset_all = zeros(numel(var_all), num_params, nruns); 
trav_wave_all = zeros(numel(var_all), num_params, nruns); 
trav_wave_all_2 = zeros(numel(var_all), num_params, nruns); 

%pattern = sprintf(...
%    'two_signal_mult_N%d_ini_state_TW_params_%s_mcsteps_%s_t_out_%s_period_%s-v%s',...
%    N, '(\d+)', '(\d+)', '(\d+)', '(\d+|Inf)', '(\d+)'); % '.' = anything
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_TW_params_%s_noise_%s_t_out_%s_period_%s-v%s',...
    N, '(\d+)', '(\d+p\d+)', '(\d+)', '(\d+|Inf)', '(\d+)'); % '.' = anything

for i=1:numel(names)
    if isempty(regexp(names{i}, pattern, 'once')) % only load files matching a certain pattern
        continue
    else
        disp(names{i});
        load( fullfile( folder, strcat(names{i}, '.mat')), 'cells_hist',...
            'save_consts_struct', 'distances', 'period', 't_out', 't_onset');
        
        [tokens, ~] = regexp(names{i}, pattern, 'tokens', 'match');
        
        %mcsteps = save_consts_struct.mcsteps;
        %idx = find(mcsteps == mcsteps_all, 1);
        noise = save_consts_struct.noise;
        idx = find(noise == noise_all, 1);
        
        idx2 = str2double(tokens{1}{1});
        %disp(idx2);
        
        %{
        if ~isempty(idx)
            filecount(idx, idx2) = filecount(idx, idx2) + 1;
            idx3 = filecount(idx, idx2);
            if idx3>nruns
                % only do up to nruns simulations
                continue
            end
            %disp(idx3);
            
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
        end
        %}
    end
end
%}
%% Load analyzed data
path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_with_noise';
fname_str = 'analyzed_data_TW_propagation_network_15_nruns_100_digits_5_old_v1';
load( fullfile(path, fname_str), 'filecount', 't_out_all', 'period_all',...
    't_onset_all', 'tmax', 'num_params', 'digits', 'trav_wave_all', 'trav_wave_all_2' );

%% Load data files 
% (Part 2)
%path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\randomized lattice';
path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_with_noise';
appendix = '_part2';
subfolder = sprintf('TW_propagation_network_%d%s', network, appendix);

folder = fullfile(path, subfolder);

listing = dir(folder);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
names = {};
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
% (Part 2)
filecount_p2 = zeros(numel(noise_all_2), num_params); 
t_out_all_p2 = zeros(numel(noise_all_2), num_params, nruns); % final times
period_all_p2 = zeros(numel(noise_all_2), num_params, nruns); % periodicity test
t_onset_all_p2 = zeros(numel(noise_all_2), num_params, nruns); 
trav_wave_all_p2 = zeros(numel(noise_all_2), num_params, nruns); 
trav_wave_all_2_p2 = zeros(numel(noise_all_2), num_params, nruns); 

%pattern = sprintf(...
%    'two_signal_mult_N%d_ini_state_TW_params_%s_mcsteps_%s_t_out_%s_period_%s-v%s',...
%    N, '(\d+)', '(\d+)', '(\d+)', '(\d+|Inf)', '(\d+)'); % '.' = anything
pattern = sprintf(...
    'two_signal_mult_N%d_ini_state_TW_params_%s_noise_%s_t_out_%s_period_%s-v%s',...
    N, '(\d+)', '(\d+p\d+)', '(\d+)', '(\d+|Inf)', '(\d+)'); % '.' = anything

for i=1:numel(names)
    if isempty(regexp(names{i}, pattern, 'once')) % only load files matching a certain pattern
        continue
    else
        disp(names{i});
        load( fullfile( folder, strcat(names{i}, '.mat')), 'cells_hist',...
            'save_consts_struct', 'distances', 'period', 't_out', 't_onset');
        
        [tokens, ~] = regexp(names{i}, pattern, 'tokens', 'match');
        
        %mcsteps = save_consts_struct.mcsteps;
        %idx = find(mcsteps == mcsteps_all, 1);
        noise = save_consts_struct.noise;
        idx = find(noise == noise_all_2, 1);
        
        idx2 = str2double(tokens{1}{1});
        %disp(idx2);
        
        if ~isempty(idx)
            filecount_p2(idx, idx2) = filecount_p2(idx, idx2) + 1;
            idx3 = filecount_p2(idx, idx2);
            if idx3>nruns
                % only do up to nruns simulations
                continue
            end
            %disp(idx3);
            
            t_out_all_p2(idx, idx2, idx3) = t_out;
            period_all_p2(idx, idx2, idx3) = period;
            t_onset_all_p2(idx, idx2, idx3) = t_onset;
            
            % test for TW
            a0 = save_consts_struct.a0;
            digits = 3;
            [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                period, t_out, distances, digits);
            trav_wave_all_p2(idx, idx2, idx3) = trav_wave;
            trav_wave_all_2_p2(idx, idx2, idx3) = trav_wave_2;
        end
        %}
    end
end

%% Merge data
% process (old) data based on the set number of simulations
t_out_all = t_out_all(:, :, 1:nruns);
period_all = period_all(:, :, 1:nruns);
t_onset_all = t_onset_all(:, :, 1:nruns);
trav_wave_all = trav_wave_all(:, :, 1:nruns);
trav_wave_all_2 = trav_wave_all_2(:, :, 1:nruns);

% merge data
noise_all_new = [noise_all noise_all_2];
filecount_new = [filecount; filecount_p2];
t_out_all_new = [t_out_all; t_out_all_p2];
period_all_new = [period_all; period_all_p2];
t_onset_all_new = [t_onset_all; t_onset_all_p2];
trav_wave_all_new = [trav_wave_all; trav_wave_all_p2];
trav_wave_all_2_new = [trav_wave_all_2; trav_wave_all_2_p2];

%{
t_out_all_p2 = zeros(numel(noise_all_2), num_params, nruns); % final times
period_all_p2 = zeros(numel(noise_all_2), num_params, nruns); % periodicity test
t_onset_all_p2 = zeros(numel(noise_all_2), num_params, nruns); 
trav_wave_all_p2 = zeros(numel(noise_all_2), num_params, nruns); 
trav_wave_all_2_p2 = zeros(numel(noise_all_2), num_params, nruns); 
%}

%% Save the analyzed data
%
subfolder = sprintf('TW_propagation_network_%d', network);
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);
save_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_with_noise';
save( fullfile(save_path, strcat(fname_str, '.mat')), 'noise_all_new',...
    'filecount_new', 't_out_all_new', 'period_all_new', 't_onset_all_new', 'tmax',...
    'num_params', 'nruns', 'digits', 'trav_wave_all_new', 'trav_wave_all_2_new');
%save( fullfile(save_path, strcat(fname_str, '.mat')), 'mcsteps_all',...
%    'filecount', 't_out_all', 'period_all', 't_onset_all', 'tmax',...
%    'num_params', 'nruns', 'digits', 'trav_wave_all', 'trav_wave_all_2');
%}
%% Load the analyzed data
subfolder = sprintf('TW_propagation_network_%d', network);
fname_str = sprintf('analyzed_data_%s_nruns_%d_digits_5', subfolder, nruns);
save_path = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_with_noise';
load( fullfile(save_path, strcat(fname_str, '.mat')), 'noise_all_new',...
    'filecount_new', 't_out_all_new', 'period_all_new', 't_onset_all_new', 'tmax',...
    'num_params', 'nruns', 'digits', 'trav_wave_all_new', 'trav_wave_all_2_new');

noise_all = noise_all_new;
filecount = filecount_new;
t_out_all = t_out_all_new;
period_all = period_all_new;
t_onset_all = t_onset_all_new;
trav_wave_all = trav_wave_all_new;
trav_wave_all_2 = trav_wave_all_2_new;
%% Analyze fraction of TWs
%trav_wave_all_mean = sum(sum(trav_wave_all, 3), 2)/(num_params*nruns);
trav_wave_all_2_mean = sum(sum(trav_wave_all_2, 3), 2)/(num_params*nruns);

h = figure;
%x_data = [0 log10(mcsteps_all(2:end)) ]; 
%x_data = log10(mcsteps_all); 
x_data = noise_all;
y_data = trav_wave_all_2_mean;

plot( x_data(1:end), y_data(1:end), 'bo-', 'LineWidth', 1.5 );
%plot( x_data, y_data, 'bo-', 'LineWidth', 1.5 );
%plot( x_data, trav_wave_all_2_mean, 'bo-', 'LineWidth', 1.5 );
ylim([0 1]);
set(gca, 'XScale', 'log');

% labels (change per case)
%xlabel('MC steps');
xlabel('Noise $\alpha$');
ylabel('Fraction TW');

% tick labels
sel_idx = [1 4 7 10 12];
%A = sprintfc('10^{%d}', log10(x_data(sel_idx)) );
A = sprintfc('10^{%d}', log10(x_data(sel_idx)) );

%{
A = cell( numel(mcsteps_all), 1 );
A{1} = '0';
A(2:end) = sprintfc('%d', mcsteps_all(2:end) );
%}
set(gca, 'FontSize', 20, 'XTick', x_data(sel_idx), 'XTickLabels', A);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_frac_TW_all_mean'));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end
%% fraction with period 15
period_15_frac = sum( sum( period_all==15, 3 ), 2)/(num_params*nruns);
period_15_mult_frac = sum( sum( mod(period_all, 15)==0, 3 ), 2)/(num_params*nruns);

h = figure;
hold on
%x_data = [0 log10(mcsteps_all(2:end)) ]; 
%x_data = log10(mcsteps_all); 
x_data = log10(noise_all);

%plot( x_data, period_15_frac, 'bo-', 'LineWidth', 1.5 );
%plot( x_data, period_15_mult_frac, 'ro-', 'LineWidth', 1.5 );
plot( x_data(1:end), period_15_frac(1:end), 'bo-', 'LineWidth', 1.5 );
plot( x_data(1:end), period_15_mult_frac(1:end), 'ro-', 'LineWidth', 1.5 );
ylim([0 1]);
set(gca, 'FontSize', 20);
legend({'T=15', 'mod(T,15)=0'}, 'Location', 'sw');

% labels (change per case)
%xlabel('MC steps');
xlabel('Noise $\alpha$');
ylabel('Fraction period 15');

% tick labels
sel_idx = [1 4 7 10 12];
A = sprintfc('10^{%d}', x_data(sel_idx) );
%A = sprintfc('%.2f', noise_all );
%{
A = cell( numel(mcsteps_all), 1 );
A{1} = '0';
A(2:end) = sprintfc('%d', mcsteps_all(2:end) );
%}
set(gca, 'FontSize', 20, 'XTick', x_data(sel_idx), 'XTickLabels', A);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_frac_period_15'));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end
%% Analyze t_out
% (1) <t_out>
%t_out_mean = mean(t_out_all, 3);

% (2) # trajectories reaching tmax
t_max_reached = (t_out_all == tmax);
t_max_count = sum(sum(t_max_reached, 3), 2);

h2 = figure;
%x_data = [0 log10(mcsteps_all(2:end)) ]; 
%x_data = log10(mcsteps_all); 
x_data = noise_all; 

plot( x_data(1:end), t_max_count(1:end)/(num_params*nruns), 'bo-', 'LineWidth', 1.5 );
set(gca, 'FontSize', 20);
ylim([0 1]);
set(gca, 'XScale', 'log');

%xlabel('MC steps');
xlabel('Noise $\alpha$');
%xlabel('MC steps');
ylabel('Fraction reaching $t_{max}$');
title(sprintf('$t_{max}=10^{%d}$', log10(tmax)));

%imagesc(p1, p2, t_max_count/nruns);
%set(gca, 'YDir', 'Normal');
%c = colorbar;
%caxis([0 1]);
%colormap('parula');
%xlabel('$$p_1$$')
%ylabel('$$p_2$$')
%ylabel(c, 'fraction');

% tick labels
sel_idx = [1 4 7 10 12];
A = sprintfc('10^{%d}', log10(x_data(sel_idx)) );
%A = sprintfc('%.2f', x_data );
%{
A = cell( numel(mcsteps_all), 1 );
A{1} = '0';
A(2:end) = sprintfc('%d', mcsteps_all(2:end) );
%}
set(gca, 'FontSize', 20, 'XTick', x_data(sel_idx), 'XTickLabels', A);

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, strcat('analyzed_data_', subfolder,...
        sprintf('_nruns_%d_digits_%d', nruns, digits), '_frac_reaching_tmax'));
    save_figure(h2, 10, 8, fname, '.pdf', qsave);
end
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