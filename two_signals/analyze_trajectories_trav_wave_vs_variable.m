%% Analyze saved trajectories for signs of travelling waves 
% Check across a range of a variable (=noise, Hill)
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');
%% Parameters
N = 225;
a0 = 1.5;
nruns = 100;
tmax = 10000;
%noise_all = [0 0.01 0.05 0.1 0.5];
hill_all = [2 5 10 100 Inf];
K12 = 9;
digits = 3;

% variable to loop over
loopvar = hill_all; %noise_all; %
loopvar_str = 'hill'; %'noise'; 

% folder for saving figures
save_path_fig = 'H:\My Documents\Multicellular automaton\figures\two_signals\travelling_wave';

%% Get all filenames
%path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution\vs_p_ini_batch3';
%path = 'D:\Multicellularity\data\two_signals\time_evolution\vs_pini_batch3';
%parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_noise\K12_9';
parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_Hill\K12_9';

names = {}; %cell(numel(loopvar), nruns);
for i1=1:numel(loopvar)
    folder = parent_folder; %fullfile(parent_folder, subfolder);
    listing = dir(folder);
    num_files = numel(listing)-2; %first two entries are not useful
    count = 0;
    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,fname,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            count = count + 1;
            names{count} = fname;
        end
    end
end

%% filename pattern
% get file info from regexp
pat1 = '\d+|Inf';
pat2 = '(\d+p\d+|Inf)';
pat3 = '\w*';

if strcmp(loopvar_str, 'noise')
    pattern = sprintf('two_signal_mult_N%d_noise_%s_K12_%d_t_out_%s_period_%s-v%s%s',...
        N, pat2, K12, pat1, pat1, pat1, pat3);
    %pattern = sprintf('two_signal_mult_N%d_initiateI0_noise_%s_t_out_%s_period_%s_%stemp-v%s',...
    %    N, pat2, pat1, pat1, pat3, pat1); % noise
elseif strcmp(loopvar_str, 'hill')
    pattern = sprintf('two_signal_mult_N%d_hill_%s_K12_%d_t_out_%s_period_%s-v%s%s',...
        N, pat2, K12, pat1, pat1, pat1, pat3); % Hill
    %pattern = sprintf('two_signal_mult_N%d_initiateI0_hill_%s_t_out_%s_period_%s_%stemp-v%s',...
    %    N, pat2, pat1, pat1, pat3, pat1); % Hill
end
%% Load raw data
% variables to store
fname_all = cell(numel(loopvar), nruns);
filecount = zeros(numel(loopvar), 1);
t_out_all = zeros(numel(loopvar), nruns); % final times
period_all = zeros(numel(loopvar), nruns); % periodicity test
t_onset_all = zeros(numel(loopvar), nruns);
trav_wave_all = zeros(numel(loopvar), nruns);
trav_wave_2_all = zeros(numel(loopvar), nruns);

error_files = {};
for i1=1:numel(names)
    %test_fname1 = 'two_signal_mult_N225_initiateI0_noise_0p50_t_out_10000_period_Inf_tmax_reached_temp-v1';
    %test_fname2 = 'two_signal_mult_N225_initiateI0_noise_0p10_t_out_3491_period_15_temp-v1';
    fname = names{i1};
    
    [~, tokens] = regexp(fname, pattern, 'match', 'tokens');
    
    % get variable from regexp tokens, match to variable list
    this_var = str2double(strrep(tokens{1}{1}, 'p', '.')); % noise / Hill
    
    idx = find(this_var == loopvar);
    
    if isempty(idx) % only load files matching a certain pattern
        continue
    else
        disp(fname);
        load(fullfile(folder, fname));
        
        %
        filecount(idx) = filecount(idx) + 1;
        idx2 = filecount(idx);
        fname_all{idx, idx2} = fname; % store filename
        t_out_all(idx, idx2) = t_out;
        period_all(idx, idx2) = period;
        t_onset_all(idx, idx2) = t_onset;
        %trav_wave_all(idx, idx2) = trav_wave;
        if period<Inf
            [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0,...
                period, t_out, digits);
            trav_wave_all(idx, idx2) = trav_wave;
            trav_wave_2_all(idx, idx2) = trav_wave_2;
        % else: not a travelling wave, no need to update trav_wave_all
        end
        %}
    end
end
%}

%% tester
%{
folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_Hill';
fname = 'two_signal_mult_N225_initiateI0_hill_100p00_t_out_1503_period_15_temp-v1';
load(fullfile(folder, fname));
if period<Inf
    [trav_wave, trav_wave_2] = travelling_wave_test(cells_hist, a0, period, t_out);
end
%}
%{
test_str = fname; % 'two_signal_mult_N225_initiateI0_hill_5p00_t_out_9495_period_10_tmax_reached_temp-v1';
[match, tokens] = regexp(test_str, pattern, 'match', 'tokens');
this_var = str2double(strrep(tokens{1}{1}, 'p', '.')); % noise / Hill
%}
% Select data to examine
%{
%select_idx = (trav_wave_all ~= trav_wave_2_all);
select_idx = trav_wave_all==1;
fname_all(select_idx);
%}
%% Save the loaded raw data
%
fname_str = sprintf('trav_wave_occur_vs_%s_K12_%d_nruns%d_round_p_I_digits_%d',...
    loopvar_str, K12, nruns, digits);
save_path = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\analyzed_data';

save(fullfile(save_path, strcat(fname_str, '.mat') ), 'loopvar', 'loopvar_str', 'filecount',...
    't_out_all', 'period_all', 't_onset_all', 'trav_wave_all', 'trav_wave_2_all') %, save_vars);
%}

%% Load processed data
%{
fname_str = sprintf('trav_wave_occur_vs_%s_K12_%d_nruns%d_round_p_I_digits_%d',...
    loopvar_str, K12, nruns, digits);
load_path = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\analyzed_data';
load(fullfile(load_path, strcat(fname_str, '.mat')));
%}

%% Fraction travelling waves
h = figure;
hold on
%---
if strcmp(loopvar_str, 'hill')
    xdata = [loopvar(1:end-1) 10^3];
    plot(xdata, sum(trav_wave_all, 2)/nruns, 'bo--', 'LineWidth', 2);
    plot(xdata, sum(trav_wave_2_all, 2)/nruns, 'ro--', 'LineWidth', 2);
    xlabel('Hill coeff. $n$')
    xtick_str = [string(loopvar(1:end-1)), "Inf"];
%---
elseif strcmp(loopvar_str, 'noise')
    xdata = [0.001 loopvar(2:end)];
    plot(xdata, sum(trav_wave_all, 2)/nruns, 'bo--', 'LineWidth', 2);
    plot(xdata, sum(trav_wave_2_all, 2)/nruns, 'ro--', 'LineWidth', 2);
    xlabel('Noise $\alpha$')
    xtick_str = string([0 loopvar(2:end)]);
end
%---
ylabel('Fraction')
title('Fraction travelling waves');
set(gca, 'FontSize', 24);
set(gca, 'XScale', 'log', 'XTick', xdata, 'XTickLabels',...
    xtick_str);
ylim([0 1]);
legend('strict crit.', 'loose crit.');

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, sprintf('frac_trav_waves_vs_%s_K12_%d_nruns_%d_p_I_%d_digits',...
        loopvar_str, K12, nruns, digits));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end

%% Filter data selectively
% multiples of period 15
frac_period_15_mult = sum((mod(period_all, sqrt(N)) == 0), 2)/nruns;

% Estimate by hand
%frac_trav_wave_manual = [0 4 33 89 69]/100;

h = figure;
hold on
if strcmp(loopvar_str, 'hill')
    %---
    xdata = [loopvar(1:end-1) 10^3];
    xlabel('Hill coeff. $n$')
    xtick_str = [string(loopvar(1:end-1)), "Inf"];
    %---
elseif strcmp(loopvar_str, 'noise')
    xdata = [0.001 loopvar(2:end)];
    xlabel('Noise $\alpha$')
    xtick_str = string([0 loopvar(2:end)]);
end
%---
plot(xdata, frac_period_15_mult, 'go--', 'LineWidth', 2);
%plot(xdata, frac_trav_wave_manual, 'mo--', 'LineWidth', 2);
ylabel('Fraction')
title('Fraction travelling waves');
set(gca, 'FontSize', 24);
set(gca, 'XScale', 'log', 'XTick', xdata, 'XTickLabels',...
    xtick_str);
ylim([0 1]);
legend('period = multiple of 15', 'Location', 'nw');

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, sprintf(...
        'frac_trav_waves_vs_%s_method2_K12_%d_nruns_%d_p_I_%d_digits',...
        loopvar_str, K12, nruns, digits));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end

%% Plot t_onset vs Hill

% process data
t_onset_mean = Inf*ones(numel(loopvar), 1);
t_onset_std = zeros(numel(loopvar), 1);
for i=1:numel(loopvar)
    idx = find(mod(period_all(i,:), sqrt(N)) == 0);
    if ~isempty(idx)
        t_onset_mean(i) = mean( t_onset_all(i, idx) );
        t_onset_std(i) = std( t_onset_all(i, idx) );
    end
end

if strcmp(loopvar_str, 'hill')
    %------
    xdata = [loopvar(1:end-1) 10^3];
    xtick_str = [string(loopvar(1:end-1)), "Inf"];
    xlabel('Hill coeff. $n$')
    xlim([1 10^3]);
    %------
elseif strcmp(loopvar_str, 'noise')
    xdata = [0.001 loopvar(2:end)];
    xtick_str = string([0 loopvar(2:end)]);
    xlabel('Noise $\alpha$')
    xlim([0.001 1]);
    %------
end

h=figure;
hold on
errorbar(xdata, t_onset_mean, t_onset_std,...
    'go--', 'LineWidth', 2);
ylabel('$t_{onset}$');
title('Onset time of trav. wave');
set(gca, 'FontSize', 24);
set(gca, 'XScale', 'log', 'XTick', xdata, 'XTickLabels',...
    xtick_str);
%xlim([2 10^3]);
%ylim([0 10^4]);
%legend('period = multiple of 15', 'manual', 'Location', 'nw');

qsave = 1;
if qsave
    fname = fullfile(save_path_fig, sprintf(...
        't_onset_trav_waves_vs_%s_K12_%d_nruns_%d_digits_%d_errorbar',...
        loopvar_str, K12, nruns, digits));
    save_figure(h, 10, 8, fname, '.pdf', qsave);
end

% Box plot
%{
idx = (mod(period_all, sqrt(N)) == 0);
box_data = t_onset_all(idx);
group_data = repmat(1:5', 1, nruns);
group_data = group_data(idx);
boxplot(box_data, group_data, 'Labels', [string(loopvar(2:end-1)), "Inf"]);
%}